# GitHub Copilot Instructions for Pairing Library Development

## Project Context
This repository contains implementations of bilinear pairing (specifically Optimal-Ate pairing) on three different elliptic curves: KSS-16, BLS-12, and BN. The implementation targets 128-bit security level and is based on research presented at IndoCrypt 2017.

The codebase is being modernized, modularized, and optimized according to the requirements in `project_requirements.md`. GitHub Copilot should assist developers in implementing these improvements while maintaining the mathematical correctness of the cryptographic operations.

## Mathematical Background
When suggesting code for this project, consider the following cryptographic and mathematical principles:

- **Bilinear Pairings**: A bilinear pairing is a map e: G₁ × G₂ → GT where G₁, G₂ are elliptic curve groups and GT is a multiplicative group. The map satisfies e(aP, bQ) = e(P,Q)^(ab) for any integers a,b.

- **Optimal Ate Pairing**: An optimization of the Ate pairing that reduces the Miller loop length, improving computational efficiency.

- **Elliptic Curve Types**:
  - **BN Curves**: Barreto-Naehrig curves with embedding degree k=12
  - **BLS Curves**: Barreto-Lynn-Scott curves with embedding degree k=12
  - **KSS Curves**: Kachisa-Schaefer-Scott curves with embedding degree k=16

- **GMP Integration**: The code uses GMP (GNU Multiple Precision Arithmetic Library) for handling large integers required in cryptographic operations.

## Code Structure Guidelines

### Modular Design
- Each curve implementation (BN, BLS, KSS) should be in its own module
- Common operations should be extracted to shared utility modules
- Clear interface definitions should separate modules

### Memory Management
- All allocated memory must be properly freed
- Use consistent patterns for initialization and cleanup
- Avoid memory leaks in error conditions

### Error Handling
- Functions should have clear error reporting
- Use consistent return codes or error structures
- Document error conditions for each function

### Naming Conventions
- Use snake_case for functions and variables
- Use UPPER_CASE for constants
- Prefix functions with module name (e.g., `bn_init`, `bls_pair`)
- Suffix types with _t (e.g., `ec_point_t`, `fp_element_t`)

## Documentation Guidelines
- Each function should have a docstring explaining:
  - Purpose
  - Parameters
  - Return values
  - Error conditions
  - Mathematical background (where relevant)
  - Performance characteristics (where relevant)

## Examples for Code Generation

### Function Template
```c
/**
 * @brief [Brief description of function]
 * 
 * [Detailed description including mathematical background]
 * 
 * @param[out] result [Description of parameter]
 * @param[in] input1 [Description of parameter]
 * @param[in] input2 [Description of parameter]
 * 
 * @return 0 on success, error code otherwise
 */
int module_function_name(type_t* result, const type_t* input1, const type_t* input2) {
    // Parameter validation
    if (!result || !input1 || !input2) {
        return ERROR_INVALID_PARAMETER;
    }
    
    // Function implementation
    
    return SUCCESS;
}
```

### Memory Management Pattern
```c
// Initialization
int module_init(module_context_t* ctx) {
    if (!ctx) {
        return ERROR_INVALID_PARAMETER;
    }
    
    // Initialize components
    
    return SUCCESS;
}

// Cleanup
void module_clear(module_context_t* ctx) {
    if (!ctx) {
        return;
    }
    
    // Free resources
}
```

## Testing Guidance
- Each function should have corresponding unit tests
- Tests should cover:
  - Normal operation
  - Edge cases
  - Error conditions
  - Known test vectors (where available)

## Performance Considerations
- Minimize memory allocations in critical paths
- Use appropriate algorithm optimizations for each curve type
- Benchmark critical operations
- Consider constant-time implementations for security-sensitive operations

## Security Best Practices
- Avoid timing side channels in cryptographic operations
- Clear sensitive data from memory after use
- Validate all inputs
- Don't expose internal state unnecessarily
- Follow secure coding guidelines for C