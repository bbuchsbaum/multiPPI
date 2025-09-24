## Micro-DSL v3.4 Source Format Grammar

This is the HUMAN-READABLE source format for v3.4. It extends v3.3 with semantic annotations
and constrained types while maintaining clarity and completeness. A separate compilation step 
produces the compressed format.

**Budget**: Target ≤ 250 lines or ≤ 2,500 tokens for optimal LLM processing.

**1. Document Structure:**
```
@pkg package_name | description
[Type Aliases section]
[Constraint Definitions section]
[Legend section if needed]
# Package Name
[Sections with entries]
[Dependencies section]
[Meta-Footer]
```

**2. Sigils (Same as v3.3):**
```
@pkg - Package declaration
@f   - Function
@d   - Data object
@x   - Re-export from another package

S3 System:
@s3g - S3 generic (UseMethod)
@s3m - S3 method (generic.class)  
@s3c - S3 class definition

S4 System:
@s4g - S4 generic (setGeneric)
@s4m - S4 method (setMethod)
@s4c - S4 class (setClass)

S7 System:
@s7g - S7 generic (new_generic)
@s7m - S7 method
@s7c - S7 class (new_class)

R6 System:
@r6c - R6 class (R6Class)
```

**3. Type System (v3.4 Enhanced with Constraints):**
```
Scalars (default): int, dbl, chr, lgl, raw, cpl
Vectors: vec<type> or type[] 
Matrices: mat<type> or mat<type,rows,cols>
Arrays: arr<type,dims>
Lists: lst<type> or lst{field:type, ...} (structured)
Data frames: df, tbl, data.table
Factors: fct, ord
Dates: Date, POSIXct, POSIXlt

Union types: type1|type2|type3
Nullable: type? (shorthand for type|NULL)
Any type: any
Ellipsis: ... or ...:type (e.g., ...:expr for NSE)

Class types: s3:classname, s4:classname, r6:classname, s7:classname

CONSTRAINED TYPES (v3.4):
Enums: chr["opt1"|"opt2"|"opt3"]
Ranges: int[min..max], dbl[min..max]
Patterns: chr[/regex/]
Exclusions: int[1..100]&!=[13,17]
References: @ref:constraint_name
```

**4. Entry Format (v3.4 Enhanced):**
```
@sigil name (param1:type1[constraint]?=default, param2:type2, ...) | Description -> return_type | return_schema
  +tag:value +tag:value
  !cov [Class1, Class2] (for generics)
  - param1 : Additional description
    @annotation:value @annotation:value
  - param2 : (constants: "a", "b", CONST) Valid values
    @requires:condition @affects:target
  - param3 : (key_funcs: func1, func2) Related functions
    @lifecycle:init @units:measurement
  ```verbatim
  # Optional verbatim R code block
  ```
```

**5. Type Aliases Section (v3.4 Enhanced):**
```markdown
## Type Aliases:
DF = df|tbl|data.table              # Standard data frame types
V<T> = vec<T>                       # Vector shorthand
Fml = s3:formula                    # Formula objects
Gg = s3:ggplot                      # ggplot2 objects
Config = lst{method:chr, opts:lst}  # Structured config
ValidPort = int[1024..65535]       # Constrained port range
```

**Standard aliases** (use these by default):
- `DF` for data frame arguments
- `Fml` for formula arguments (not `fml`)
- `V<T>` for vectors when brevity helps

**6. Constraint Definitions (v3.4 New):**
```markdown
## Constraint Definitions:
@constraint positive_weights | Positive numeric weights
  type: vec<dbl>
  validates: all(. > 0)
  length: @env:nrow(data)
  
@constraint valid_identifier | Valid R identifier
  type: chr
  pattern: /^[a-zA-Z_][a-zA-Z0-9_.]*$/
  not_reserved: TRUE
```

**7. Class Documentation (v3.4 Enhanced):**
```
@s4c ClassName | One-line description
  - slots: name1 (type1[constraint]) @annotation:value
           name2 (type2) @lifecycle:init @immutable
  - extends: ParentClass
  - validity: Description of validity rules
  
@r6c ClassName | One-line description  
  - fields: field1 (type1[constraint]) @purpose:role
            field2 (type2) @lazy @cached
  - methods: method1 (args) -> ret_type
             method2 (args) -> ret_type  
  - inherits: ParentClass
```

**8. Metadata Tags (v3.3 + v3.4 additions):**
```
+family:group_name           # Function family
+pipe:in|out                # Pipe compatibility (in, out, or both)
+nse:param1,param2          # Parameters using NSE
+side:effect[details]       # Side effects with sub-facets
  - fs[read|write|delete]   # File system operations
  - plot[device|file]       # Graphics output
  - console[print|message|warning]  # Console output
  - network[http|socket|download]   # Network operations
  - options[get|set|env]    # Global options/environment
  - db[read|write|query]    # Database operations
+perf:O(complexity)         # Performance complexity
+mem:usage                  # Memory usage pattern
+compute:intensity          # Computational intensity
+deprecated:replacement     # Deprecation with suggested alternative
+wraps:function            # This function wraps another
+calls:func1,func2         # Functions called internally
+see:related1,related2     # Related functions to consider
+parallel:capable          # Can use parallel processing (v3.4)
+deterministic:false       # Non-deterministic results (v3.4)
+pure:false               # Has side effects (v3.4)
```

**9. Semantic Annotations (v3.4 New):**
```
BEHAVIORAL:
@controls:aspect          # Parameter controls specific behavior
@affects:target          # Changes affect another component  
@modifies:target         # Directly modifies target

DEPENDENCY:
@requires:condition      # Prerequisite condition
@conflicts:parameter     # Mutually exclusive with
@extends:base           # Extends functionality

VALIDATION:
@validates:constraint    # Validation rule
@range:[min,max]        # Numeric range
@length:constraint      # Length requirement
@pattern:regex          # Pattern matching

SEMANTIC ROLE:
@purpose:role           # Semantic purpose
@units:measurement      # Physical/logical units
@example:value          # Example values
@default-reason:why     # Why this default

LIFECYCLE:
@lifecycle:stage        # When relevant (init|config|runtime|cleanup)
@immutable             # Cannot be modified
@cached                # Result is cached
@lazy                  # Evaluated on demand

CONDITIONAL:
@when:condition        # Conditional applicability
@implies:consequence   # Logical implication
@if:cond @then:result  # If-then constraints
```

**10. Structured Return Types (v3.4 New):**
```
-> lst{
  field1: type1 @annotation,
  field2: type2[constraint] @annotation,
  nested: lst{
    subfield: type
  }
}
```

**11. Example Entry (v3.4):**
```markdown
@f analyze_model (
  model:s3:lm,
  type:chr["summary"|"anova"|"diagnostics"]?="summary",
  conf.level:dbl[0.5..0.99]?=0.95
) | Analyze fitted model -> lst{
  statistics: df @purpose:results,
  plots: lst<s3:ggplot>? @when:type="diagnostics",
  interpretation: chr @purpose:summary
}
  +family:analysis +compute:light
  - model : @requires:fitted @validates:has-residuals
  - type : @controls:output-format @affects:return-structure
  - conf.level : @purpose:confidence @affects:statistics.ci
```

**12. Conditional Constraints (v3.4):**
```markdown
@f process_data (
  data:df,
  method:chr["scale"|"center"|"none"]?="none",
  scale.center:lgl?=TRUE,
  scale.scale:lgl?=TRUE
) | Process data with scaling options -> df
  - method : @controls:processing
  - scale.center : @when:method="scale" @requires:TRUE
                   @when:method="center" @implies:scale.scale=FALSE
  - scale.scale : @when:method="scale" @default:TRUE
                  @conflicts:method="center"
```

**13. Best Practices (v3.4):**
- Use specific sigils (@s3g not @g)
- Always specify vector vs scalar types
- Use standard type aliases (DF, Fml, V<T>)
- Add constraints from match.arg/stopifnot/checks
- Keep !cov lists short (3-6 classes max)
- Document semantic relationships concisely
- Use structured types for complex returns
- Define reusable constraints with @constraint
- Include conditional logic with @when/@implies
- Group related functions with +family tags
- Mark side effects with detailed sub-facets
- Stay within budget (≤250 lines)

**14. Meta-Footer:**
```markdown
---
## Meta-Footer
- Micro-DSL Version: v3.4-source
- Package: {pkg} (Version: X.Y.Z)
- Generated: [ISO-8601 timestamp]
- Features: types[constrained] sigils[specific] metadata[rich] semantics[annotated]
- Coverage: {n_documented_exports} / {n_total_exports} exports
- Provenance: exports[NAMESPACE], enums[match.arg/switch], constraints[assertions/checks]
```

**15. Export Detection Priority:**
1. NAMESPACE file: `export()`, `S3method()`, `exportClasses()`, `exportMethods()`, `exportPattern()`
2. Roxygen tags: `@export` in documentation
3. If neither present: skip the symbol (do not guess or include)

**16. Inference Heuristics (apply silently):**
- Type from defaults: TRUE/FALSE → lgl, "text" → chr, 1L → int, 1.0 → dbl
- Common patterns: data/df/tbl → DF, formula → Fml, weights → vec<dbl>
- Enums: match.arg(x, c("a","b")) → chr["a"|"b"]
- Ranges: stopifnot(x >= 0 && x <= 1) → dbl[0..1]
- Side effects: file.* → fs, plot/ggplot → plot, message/cat → console
- Determinism: runif/sample/rnorm → +deterministic:false

---

```markdown
@pkg fmrihrf | Hemodynamic Response Functions for fMRI Data Analysis

## Type Aliases:
DF = df|tbl|data.table
V<T> = vec<T>
Fml = s3:formula
Gg = s3:ggplot

## Constraint Definitions:
@constraint positive_weights | Positive numeric weights
  type: vec<dbl>
  validates: all(. > 0)

# fmrihrf

## 1. Core

@f HRF (fun:chr|function, name:chr, nbasis:int?=1, span:dbl?=24, param_names:chr[]?=NULL) | Create HRF object -> s3:HRF
  +family:hrf +pipe:in|out
  - fun : @purpose:hrf-function
  - name : @purpose:hrf-name
  - nbasis : @purpose:number-of-basis-functions
  - span : @purpose:hrf-span
  - param_names : @purpose:parameter-names

@f HRF_BSPLINE (nbasis:int?=5, span:dbl?=24) | B-spline HRF -> s3:HRF
  +family:hrf +pipe:in|out
  - nbasis : @purpose:number-of-basis-functions
  - span : @purpose:hrf-span

@f HRF_FIR (nbasis:int?=12, span:dbl?=24) | FIR HRF -> s3:HRF
  +family:hrf +pipe:in|out
  - nbasis : @purpose:number-of-basis-functions
  - span : @purpose:hrf-span

@f HRF_GAMMA (t:vec<dbl>, shape:dbl?=6, rate:dbl?=1) | Gamma HRF -> vec<dbl>
  +family:hrf +pipe:in|out
  - t : @purpose:time-points
  - shape : @purpose:shape-parameter
  - rate : @purpose:rate-parameter

@f HRF_GAUSSIAN (t:vec<dbl>, mean:dbl?=6, sd:dbl?=2) | Gaussian HRF -> vec<dbl>
  +family:hrf +pipe:in|out
  - t : @purpose:time-points
  - mean : @purpose:mean-parameter
  - sd : @purpose:sd-parameter

@f HRF_SPMG1 (t:vec<dbl>, P1:dbl?=5, P2:dbl?=15, A1:dbl?=0.0833) | SPM canonical HRF -> vec<dbl>
  +family:hrf +pipe:in|out
  - t : @purpose:time-points
  - P1 : @purpose:shape-parameter
  - P2 : @purpose:shape-parameter
  - A1 : @purpose:amplitude-parameter

@f HRF_SPMG2 (t:vec<dbl>, P1:dbl?=5, P2:dbl?=15, A1:dbl?=0.0833) | SPM canonical HRF with temporal derivative -> mat<dbl>
  +family:hrf +pipe:in|out
  - t : @purpose:time-points
  - P1 : @purpose:shape-parameter
  - P2 : @purpose:shape-parameter
  - A1 : @purpose:amplitude-parameter

@f HRF_SPMG3 (t:vec<dbl>, P1:dbl?=5, P2:dbl?=15, A1:dbl?=0.0833) | SPM canonical HRF with temporal and dispersion derivatives -> mat<dbl>
  +family:hrf +pipe:in|out
  - t : @purpose:time-points
  - P1 : @purpose:shape-parameter
  - P2 : @purpose:shape-parameter
  - A1 : @purpose:amplitude-parameter

@f acquisition_onsets (x:s3:sampling_frame) | Get fMRI acquisition onset times -> vec<dbl>
  +family:sampling +pipe:in|out
  - x : @purpose:sampling-frame

@f amplitudes (x:s3:Reg) | Get amplitudes from an object -> vec<dbl>
  +family:regressor +pipe:in|out
  - x : @purpose:regressor-object

@f as_hrf (f:function, name:chr, nbasis:int?=1, span:dbl?=24, params:lst?=list()) | Turn any function into an HRF object -> s3:HRF
  +family:hrf +pipe:in|out
  - f : @purpose:hrf-function
  - name : @purpose:hrf-name
  - nbasis : @purpose:number-of-basis-functions
  - span : @purpose:hrf-span
  - params : @purpose:hrf-parameters

@f bind_basis (...:s3:HRF) | Bind HRFs into a basis set -> s3:HRF
  +family:hrf +pipe:in|out
  - ... : @purpose:hrf-objects

@f block_hrf (hrf:s3:HRF, width:dbl, precision:dbl?=0.1, half_life:dbl?=Inf, summate:lgl?=TRUE, normalize:lgl?=FALSE) | Create a blocked HRF object -> s3:HRF
  +family:hrf +pipe:in|out
  - hrf : @purpose:hrf-object
  - width : @purpose:block-width
  - precision : @purpose:sampling-precision
  - half_life : @purpose:exponential-decay
  - summate : @purpose:summate-responses
  - normalize : @purpose:normalize-response

@f deriv (x:s3:HRF, t:vec<dbl>, ...) | Compute derivatives of HRF functions -> vec<dbl>|mat<dbl>
  +family:hrf +pipe:in|out
  - x : @purpose:hrf-object
  - t : @purpose:time-points

@f durations (x:s3:Reg) | Get durations of an object -> vec<dbl>
  +family:regressor +pipe:in|out
  - x : @purpose:regressor-object

@f empirical_hrf (t:vec<dbl>, y:vec<dbl>, name:chr?="empirical_hrf") | Generate an empirical HRF -> s3:HRF
  +family:hrf +pipe:in|out
  - t : @purpose:time-points
  - y : @purpose:hrf-values
  - name : @purpose:hrf-name

@f evaluate (x:s3:HRF|s3:Reg, grid:vec<dbl>, ...) | Evaluate an HRF or regressor object -> vec<dbl>|mat<dbl>
  +family:evaluation +pipe:in|out
  - x : @purpose:hrf-or-regressor
  - grid : @purpose:time-points

@f gen_hrf (hrf:chr|function|s3:HRF, lag:dbl?=0, width:dbl?=0, precision:dbl?=0.1, half_life:dbl?=Inf, summate:lgl?=TRUE, normalize:lgl?=FALSE, name:chr?=NULL, span:dbl?=NULL, ...) | Construct an HRF instance using decorators -> s3:HRF
  +family:hrf +pipe:in|out
  - hrf : @purpose:hrf-object
  - lag : @purpose:temporal-lag
  - width : @purpose:block-width
  - precision : @purpose:sampling-precision
  - half_life : @purpose:exponential-decay
  - summate : @purpose:summate-responses
  - normalize : @purpose:normalize-response
  - name : @purpose:hrf-name
  - span : @purpose:hrf-span

@f getHRF (name:chr, nbasis:int?=5, span:dbl?=24, lag:dbl?=0, width:dbl?=0, summate:lgl?=TRUE, normalize:lgl?=FALSE, ...) | Get HRF by name -> s3:HRF
  +family:hrf +pipe:in|out
  - name : @purpose:hrf-name
  - nbasis : @purpose:number-of-basis-functions
  - span : @purpose:hrf-span
  - lag : @purpose:temporal-lag
  - width : @purpose:block-width
  - summate : @purpose:summate-responses
  - normalize : @purpose:normalize-response

@f hrf_from_coefficients (hrf:s3:HRF, h:vec<dbl>, name:chr?=NULL, ...) | Combine HRF basis with coefficients -> s3:HRF
  +family:hrf +pipe:in|out
  - hrf : @purpose:hrf-object
  - h : @purpose:coefficients
  - name : @purpose:hrf-name

@f lag_hrf (hrf:s3:HRF, lag:dbl) | Lag an HRF object -> s3:HRF
  +family:hrf +pipe:in|out
  - hrf : @purpose:hrf-object
  - lag : @purpose:temporal-lag

@f list_available_hrfs (details:lgl?=FALSE) | List all available HRFs -> df
  +family:hrf +pipe:in|out
  - details : @purpose:include-details

@f make_hrf (basis:chr|function|s3:HRF, lag:dbl, nbasis:int?=1) | Create an HRF from a basis specification -> s3:HRF
  +family:hrf +pipe:in|out
  - basis : @purpose:hrf-basis
  - lag : @purpose:temporal-lag
  - nbasis : @purpose:number-of-basis-functions

@f nbasis (x:s3:HRF|s3:Reg, ...) | Number of basis functions -> int
  +family:hrf +pipe:in|out
  - x : @purpose:hrf-or-regressor

@f neural_input (x:s3:Reg, start:dbl?=0, end:dbl?=NULL, resolution:dbl?=0.33, ...) | Generate neural input function from event timing -> lst{time:vec<dbl>, neural_input:vec<dbl>}
  +family:regressor +pipe:in|out
  - x : @purpose:regressor-object
  - start : @purpose:start-time
  - end : @purpose:end-time
  - resolution : @purpose:temporal-resolution

@f normalise_hrf (hrf:s3:HRF) | Normalise an HRF object -> s3:HRF
  +family:hrf +pipe:in|out
  - hrf : @purpose:hrf-object

@f onsets (x:s3:Reg) | Get event onsets from an object -> vec<dbl>
  +family:regressor +pipe:in|out
  - x : @purpose:regressor-object

@f penalty_matrix (x:s3:HRF, order:int?=2, ...) | Generate penalty matrix for regularization -> mat<dbl>
  +family:hrf +pipe:in|out
  - x : @purpose:hrf-object
  - order : @purpose:penalty-order

@f plot (x:s3:HRF, ...) | Plot an HRF object -> invisible<any>
  +family:hrf +pipe:in|out +side:plot[device]
  - x : @purpose:hrf-object

@f reconstruction_matrix (hrf:s3:HRF, sframe:s3:sampling_frame|vec<dbl>, ...) | Reconstruction matrix for an HRF basis -> mat<dbl>
  +family:hrf +pipe:in|out
  - hrf : @purpose:hrf-object
  - sframe : @purpose:sampling-frame-or-time-points

@f regressor (onsets:vec<dbl>, hrf:s3:HRF?=HRF_SPMG1, duration:vec<dbl>?=0, amplitude:vec<dbl>?=1, span:dbl?=40, summate:lgl?=TRUE) | Construct a regressor object -> s3:Reg
  +family:regressor +pipe:in|out
  - onsets : @purpose:event-onsets
  - hrf : @purpose:hrf-object
  - duration : @purpose:event-durations
  - amplitude : @purpose:event-amplitudes
  - span : @purpose:hrf-span
  - summate : @purpose:summate-responses

@f regressor_design (onsets:vec<dbl>, fac:fct, block:int[], sframe:s3:sampling_frame, hrf:s3:HRF?=HRF_SPMG1, duration:vec<dbl>?=0, amplitude:vec<dbl>?=1, span:dbl?=40, precision:dbl?=0.33, method:chr["conv"|"fft"|"Rconv"|"loop"]?="conv", sparse:lgl?=FALSE, summate:lgl?=TRUE) | Build a design matrix from block-wise onsets -> mat<dbl>|sparse
  +family:regressor +pipe:in|out
  - onsets : @purpose:event-onsets
  - fac : @purpose:condition-factor
  - block : @purpose:block-identifiers
  - sframe : @purpose:sampling-frame
  - hrf : @purpose:hrf-object
  - duration : @purpose:event-durations
  - amplitude : @purpose:event-amplitudes
  - span : @purpose:hrf-span
  - precision : @purpose:sampling-precision
  - method : @purpose:evaluation-method
  - sparse : @purpose:return-sparse-matrix
  - summate : @purpose:summate-responses

@f regressor_set (onsets:vec<dbl>, fac:fct, hrf:s3:HRF?=HRF_SPMG1, duration:vec<dbl>?=0, amplitude:vec<dbl>?=1, span:dbl?=40, summate:lgl?=TRUE) | Construct a regressor set -> s3:RegSet
  +family:regressor +pipe:in|out
  - onsets : @purpose:event-onsets
  - fac : @purpose:condition-factor
  - hrf : @purpose:hrf-object
  - duration : @purpose:event-durations
  - amplitude : @purpose:event-amplitudes
  - span : @purpose:hrf-span
  - summate : @purpose:summate-responses

@f sampling_frame (blocklens:vec<int>, TR:vec<dbl>, start_time:vec<dbl>?=TR/2, precision:dbl?=0.1) | Create a sampling frame -> s3:sampling_frame
  +family:sampling +pipe:in|out
  - blocklens : @purpose:block-lengths
  - TR : @purpose:repetition-time
  - start_time : @purpose:start-time
  - precision : @purpose:sampling-precision

@f shift (x:s3:Reg, shift_amount:dbl, ...) | Shift a time series object -> s3:Reg
  +family:regressor +pipe:in|out
  - x : @purpose:regressor-object
  - shift_amount : @purpose:shift-amount

@f single_trial_regressor (onsets:dbl, hrf:s3:HRF?=HRF_SPMG1, duration:dbl?=0, amplitude:dbl?=1, span:dbl?=24) | Create a single trial regressor -> s3:Reg
  +family:regressor +pipe:in|out
  - onsets : @purpose:event-onset
  - hrf : @purpose:hrf-object
  - duration : @purpose:event-duration
  - amplitude : @purpose:event-amplitude
  - span : @purpose:hrf-span

## Dependencies
- Imports: Rcpp, assertthat, purrr, stats, Matrix, cli, memoise, numDeriv, splines, pracma
- Suggests: testthat, knitr, rmarkdown, ggplot2, dplyr, tidyr, viridis, scales, microbenchmark

---
## Meta-Footer
- Micro-DSL Version: v3.4-source
- Package: fmrihrf (Version: 0.1.0)
- Generated: [ISO-8601]
- Features: types[constrained] sigils[specific] metadata[rich] semantics[annotated]
- Coverage: 50 / 50 exports
- Provenance: exports[NAMESPACE], enums[match.arg/switch], constraints[assertions/checks]
```