#include <immintrin.h> // For AVX-512 intrinsics
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Define constants for the number of 64-bit limbs
#define NUM_LIMBS_512BIT 8   // 512 bits = 8 * 64 bits
#define NUM_LIMBS_1024BIT 16 // 1024 bits = 16 * 64 bits

// --- Helper Functions for String to Limbs and Limbs to String ---

// Function to convert a hex string to an array of 64-bit limbs
// The limbs are stored in little-endian order (least significant limb at index
// 0)
void hex_to_limbs(const char *hex_str, uint64_t *limbs, int num_limbs) {
  memset(limbs, 0, num_limbs * sizeof(uint64_t));
  int hex_len = strlen(hex_str);
  int limb_idx = 0;
  int hex_chars_per_limb = 16; // 64 bits = 16 hex characters

  for (int i = hex_len - 1; i >= 0; i -= hex_chars_per_limb) {
    if (limb_idx >= num_limbs)
      break;

    char temp_hex[hex_chars_per_limb + 1];
    int start_idx = i - hex_chars_per_limb + 1;
    if (start_idx < 0)
      start_idx = 0;
    int len = i - start_idx + 1;

    strncpy(temp_hex, hex_str + start_idx, len);
    temp_hex[len] = '\0';

    limbs[limb_idx++] = strtoull(temp_hex, NULL, 16);
  }
}

// Function to convert an array of 64-bit limbs to a hex string
// Assumes limbs are in little-endian order.
// Returns a dynamically allocated string, caller must free.
char *limbs_to_hex(const uint64_t *limbs, int num_limbs) {
  // Determine the exact length needed (leading zeros trimmed for the whole
  // number)
  int actual_limbs = num_limbs;
  while (actual_limbs > 0 && limbs[actual_limbs - 1] == 0) {
    actual_limbs--;
  }

  if (actual_limbs == 0) {
    return strdup("0");
  }

  // Each limb is 16 hex chars + 1 for null terminator per limb for safety.
  // Plus 1 for '0x' if needed, but we'll omit for raw hex.
  char *hex_str = (char *)malloc(actual_limbs * 16 + 1);
  if (!hex_str) {
    perror("malloc failed");
    return NULL;
  }
  hex_str[0] = '\0'; // Initialize as empty string

  char temp_limb_str[17]; // 16 chars + null terminator

  for (int i = actual_limbs - 1; i >= 0; i--) {
    sprintf(temp_limb_str, "%016llx", (unsigned long long)limbs[i]);
    strcat(hex_str, temp_limb_str);
  }

  // Trim leading zeros from the entire string if present (except if it's "0")
  if (hex_str[0] == '0' && strlen(hex_str) > 1) {
    char *non_zero_start = hex_str;
    while (*non_zero_start == '0' && *(non_zero_start + 1) != '\0') {
      non_zero_start++;
    }
    if (non_zero_start != hex_str) {
      memmove(hex_str, non_zero_start, strlen(non_zero_start) + 1);
    }
  }

  return hex_str;
}

// --- AVX-512 512-bit Integer Multiplication ---

// Function to multiply two 512-bit integers (represented by 8 uint64_t limbs)
// The result is a 1024-bit integer (16 uint64_t limbs).
// All limb arrays are little-endian.
void multiply_512bit_avx512(const uint64_t *a_limbs, const uint64_t *b_limbs,
                            uint64_t *result_limbs) {
  // Initialize result to zero
  memset(result_limbs, 0, NUM_LIMBS_1024BIT * sizeof(uint64_t));

  // Load one operand (a_limbs) into a temporary array for easier indexing with
  // AVX-512 We'll process 8 products at a time. The main loop will iterate
  // through the 'a' limbs. For each a_limb, we multiply it by all 8 b_limbs.
  // Then we add these partial products to the correct positions in
  // result_limbs.

  // Outer loop: Iterate through each limb of 'a'
  for (int i = 0; i < NUM_LIMBS_512BIT; ++i) {
    uint64_t current_a_limb = a_limbs[i];

    // If current_a_limb is zero, skip this iteration as it won't contribute
    if (current_a_limb == 0) {
      continue;
    }

    // Broadcast current_a_limb to all elements of a ZMM register
    __m512i va = _mm512_set1_epi64(current_a_limb);

    // Inner loop: Multiply va with batches of b_limbs
    // We'll process b_limbs in chunks of 8.
    // For 8 b_limbs, we need one __m512i load.

    // Load 8 b_limbs into a vector register
    __m512i vb = _mm512_loadu_epi64(b_limbs);

    // Perform parallel 64x64-bit multiplication
    // vp1_lo holds the lower 64 bits of each product
    __m512i vp_lo = _mm512_mullo_epi64(va, vb);
    // vp1_hi holds the upper 64 bits of each product (requires AVX512DQ)
    __m512i vp_hi = _mm512_mulhi_epi64(va, vb);

    // Store partial products to an intermediate buffer
    uint64_t temp_prod_lo[NUM_LIMBS_512BIT];
    uint64_t temp_prod_hi[NUM_LIMBS_512BIT];
    _mm512_storeu_epi64(temp_prod_lo, vp_lo);
    _mm512_storeu_epi64(temp_prod_hi, vp_hi);

    uint64_t carry = 0;

    // Add partial products to the result_limbs, handling carries
    // This is the "grade-school" addition with carry propagation
    for (int j = 0; j < NUM_LIMBS_512BIT; ++j) {
      uint64_t sum_low = result_limbs[i + j] + temp_prod_lo[j] + carry;
      carry = (sum_low < result_limbs[i + j] || sum_low < temp_prod_lo[j] ||
               (sum_low == result_limbs[i + j] && temp_prod_lo[j] == 0 &&
                carry == 1))
                  ? 1
                  : 0; // Check for overflow
      result_limbs[i + j] = sum_low;

      carry += temp_prod_hi[j]; // Add high part of product to carry
    }

    // Propagate the final carry from this row of multiplication
    int k = i + NUM_LIMBS_512BIT;
    while (carry != 0) {
      uint64_t sum_carry = result_limbs[k] + carry;
      carry = (sum_carry < result_limbs[k]) ? 1 : 0; // Check for overflow
      result_limbs[k] = sum_carry;
      k++;
      if (k >= NUM_LIMBS_1024BIT) {
        // This shouldn't happen for 512*512=1024 bit product
        // unless there's an error in logic or number of limbs is too small
        fprintf(stderr,
                "Error: Result buffer overflow during carry propagation.\n");
        break;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <hex_string_num1> <hex_string_num2>\n", argv[0]);
    return 1;
  }

  const char *s_num1 = argv[1];
  const char *s_num2 = argv[2];

  uint64_t a_limbs[NUM_LIMBS_512BIT];
  uint64_t b_limbs[NUM_LIMBS_512BIT];
  uint64_t result_limbs[NUM_LIMBS_1024BIT];

  printf("Input Num1 (hex): %s\n", s_num1);
  printf("Input Num2 (hex): %s\n", s_num2);

  hex_to_limbs(s_num1, a_limbs, NUM_LIMBS_512BIT);
  hex_to_limbs(s_num2, b_limbs, NUM_LIMBS_512BIT);

  printf("\nPerforming 512-bit multiplication using AVX-512 intrinsics...\n");
  multiply_512bit_avx512(a_limbs, b_limbs, result_limbs);

  char *result_str = limbs_to_hex(result_limbs, NUM_LIMBS_1024BIT);
  if (result_str) {
    printf("Result (hex): %s\n", result_str);
    free(result_str);
  } else {
    fprintf(stderr, "Failed to convert result to hex string.\n");
    return 1;
  }

  return 0;
}