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

```
@pkg fmridesign | Design Matrix Construction for fMRI Analysis

## Type Aliases:
DF = df|tbl|data.table
V<T> = vec<T>
Fml = s3:formula
Gg = s3:ggplot

## Constraint Definitions:
@constraint positive_weights | Positive numeric weights
  type: vec<dbl>
  validates: all(. > 0)

# fmridesign

## 1. Core

@f BSpline (x:vec<dbl>, degree:int) | B-spline basis -> s3:BSpline
  +family:basis +pipe:in +nse:x
  - x : @purpose:input
  - degree : @purpose:degree

@f Fcontrasts (x:s3:event_model, ...) | Compute F-contrasts -> lst<mat<dbl>>
  +family:contrast +pipe:in
  - x : @purpose:model

@f HRF (basis:chr["spmg1"|"spmg2"|"spmg3"|"gamma"|"gaussian"|"bspline"|"fir"], lag:dbl?=0, nbasis:int?=1) | Create HRF -> s3:HRF
  +family:hrf +pipe:in
  - basis : @purpose:basis
  - lag : @purpose:lag
  - nbasis : @purpose:number-of-basis

@f Ident (...:any) | Identity basis -> s3:Ident
  +family:basis +pipe:in +nse:...
  - ... : @purpose:variables

@f Poly (x:vec<dbl>, degree:int) | Polynomial basis -> s3:Poly
  +family:basis +pipe:in +nse:x
  - x : @purpose:input
  - degree : @purpose:degree

@f RobustScale (x:vec<dbl>) | Robust scaling -> s3:RobustScale
  +family:basis +pipe:in +nse:x
  - x : @purpose:input

@f Scale (x:vec<dbl>) | Z-score scaling -> s3:Scale
  +family:basis +pipe:in +nse:x
  - x : @purpose:input

@f ScaleWithin (x:vec<dbl>, g:chr) | Z-score within groups -> s3:ScaleWithin
  +family:basis +pipe:in +nse:x,g
  - x : @purpose:input
  - g : @purpose:group

@f Standardized (x:vec<dbl>) | Standardize -> s3:Standardized
  +family:basis +pipe:in +nse:x
  - x : @purpose:input

@f baseline (degree:int?=1, basis:chr["constant"|"poly"|"bs"|"ns"]?="constant", name:chr?=NULL, intercept:chr["runwise"|"global"|"none"]?="runwise") | Create baseline spec -> s3:baselinespec
  +family:baseline +pipe:in
  - degree : @purpose:degree
  - basis : @purpose:basis
  - name : @purpose:name
  - intercept : @purpose:intercept

@f baseline_model (basis:chr["constant"|"poly"|"bs"|"ns"]?="constant", degree:int?=1, sframe:s3:sampling_frame, intercept:chr["runwise"|"global"|"none"]?="runwise", nuisance_list:lst<mat<dbl>>?=NULL) | Construct baseline model -> s3:baseline_model
  +family:baseline +pipe:in
  - basis : @purpose:basis
  - degree : @purpose:degree
  - sframe : @purpose:sampling-frame
  - intercept : @purpose:intercept
  - nuisance_list : @purpose:nuisance

@f block (x:any) | Create block variable -> s3:blockspec
  +family:block +pipe:in
  - x : @purpose:input

@f column_contrast (pattern_A:chr, pattern_B:chr?=NULL, name:chr, where:Fml?=NULL) | Define column contrast -> s3:column_contrast_spec
  +family:contrast +pipe:in
  - pattern_A : @purpose:pattern
  - pattern_B : @purpose:pattern
  - name : @purpose:name
  - where : @purpose:subset

@f contrast (form:Fml, name:chr, where:Fml?=NULL) | Define linear contrast -> s3:contrast_formula_spec
  +family:contrast +pipe:in
  - form : @purpose:formula
  - name : @purpose:name
  - where : @purpose:subset

@f contrast_set (...:s3:contrast_spec) | Create set of contrasts -> s3:contrast_set
  +family:contrast +pipe:in +nse:...
  - ... : @purpose:contrasts

@f covariate (...:any, data:DF, id:chr?=NULL, prefix:chr?=NULL, subset:Fml?=NULL) | Construct covariate term -> s3:covariatespec
  +family:covariate +pipe:in +nse:...
  - ... : @purpose:variables
  - data : @purpose:data
  - id : @purpose:id
  - prefix : @purpose:prefix
  - subset : @purpose:subset

@f event_model (formula_or_list:Fml|lst<s3:hrfspec>, data:DF, block:Fml|vec<int>, sampling_frame:s3:sampling_frame, durations:vec<dbl>?=0, drop_empty:lgl?=TRUE, precision:dbl?=0.3, parallel:lgl?=FALSE, progress:lgl?=FALSE) | Construct event model -> s3:event_model
  +family:model +pipe:in
  - formula_or_list : @purpose:specification
  - data : @purpose:data
  - block : @purpose:block
  - sampling_frame : @purpose:sampling-frame
  - durations : @purpose:durations
  - drop_empty : @purpose:drop-empty
  - precision : @purpose:precision
  - parallel : @purpose:parallel
  - progress : @purpose:progress

@f hrf (...:any, basis:chr|s3:HRF?="spmg1", onsets:vec<dbl>?=NULL, durations:vec<dbl>?=NULL, prefix:chr?=NULL, subset:Fml?=NULL, precision:dbl?=0.3, nbasis:int?=1, contrasts:s3:contrast_spec|s3:contrast_set?=NULL, id:chr?=NULL, name:chr?=NULL, lag:dbl?=0, summate:lgl?=TRUE) | Hemodynamic regressor spec -> s3:hrfspec
  +family:hrf +pipe:in +nse:...
  - ... : @purpose:variables
  - basis : @purpose:basis
  - onsets : @purpose:onsets
  - durations : @purpose:durations
  - prefix : @purpose:prefix
  - subset : @purpose:subset
  - precision : @purpose:precision
  - nbasis : @purpose:nbasis
  - contrasts : @purpose:contrasts
  - id : @purpose:id
  - name : @purpose:name
  - lag : @purpose:lag
  - summate : @purpose:summate

@f interaction_contrast (A:Fml, name:chr, where:Fml?=NULL) | Define interaction contrast -> s3:interaction_contrast_spec
  +family:contrast +pipe:in
  - A : @purpose:formula
  - name : @purpose:name
  - where : @purpose:subset

@f one_against_all_contrast (levels:vec<chr>, facname:chr, where:Fml?=NULL) | Define one-against-all contrast -> s3:contrast_set
  +family:contrast +pipe:in
  - levels : @purpose:levels
  - facname : @purpose:factor
  - where : @purpose:subset

@f oneway_contrast (A:Fml, name:chr, where:Fml?=NULL) | Define one-way contrast -> s3:oneway_contrast_spec
  +family:contrast +pipe:in
  - A : @purpose:formula
  - name : @purpose:name
  - where : @purpose:subset

@f pair_contrast (A:Fml, B:Fml, name:chr, where:Fml?=NULL) | Define pair contrast -> s3:pair_contrast_spec
  +family:contrast +pipe:in
  - A : @purpose:formula
  - B : @purpose:formula
  - name : @purpose:name
  - where : @purpose:subset

@f poly_contrast (A:Fml, name:chr, where:Fml?=NULL, degree:int?=1, value_map:lst<dbl>?=NULL) | Define polynomial contrast -> s3:poly_contrast_spec
  +family:contrast +pipe:in
  - A : @purpose:formula
  - name : @purpose:name
  - where : @purpose:subset
  - degree : @purpose:degree
  - value_map : @purpose:value-map

@f sliding_window_contrasts (levels:vec<chr>, facname:chr, window_size:int?=2, where:Fml?=NULL, name_prefix:chr?="win") | Define sliding-window contrasts -> s3:contrast_set
  +family:contrast +pipe:in
  - levels : @purpose:levels
  - facname : @purpose:factor
  - window_size : @purpose:window-size
  - where : @purpose:subset
  - name_prefix : @purpose:prefix

@f unit_contrast (A:Fml, name:chr, where:Fml?=NULL) | Define unit contrast -> s3:unit_contrast_spec
  +family:contrast +pipe:in
  - A : @purpose:formula
  - name : @purpose:name
  - where : @purpose:subset

@f validate_contrasts (x:s3:event_model|mat<dbl>, weights:mat<dbl>|lst<mat<dbl>>?=NULL, tol:dbl?=1e-8) | Validate contrast weights -> df
  +family:contrast +pipe:in
  - x : @purpose:model
  - weights : @purpose:weights
  - tol : @purpose:tolerance

## 2. S3 Methods

@s3m Fcontrasts.convolved_term (x:s3:convolved_term, ...) | Compute F-contrasts for convolved term -> lst<mat<dbl>>
  +family:contrast +pipe:in
  !cov [convolved_term]

@s3m Fcontrasts.event_model (x:s3:event_model, ...) | Compute F-contrasts for event model -> lst<mat<dbl>>
  +family:contrast +pipe:in
  !cov [event_model]

@s3m contrast_weights.contrast_formula_spec (x:s3:contrast_formula_spec, term:s3:event_term, ...) | Compute contrast weights -> lst
  +family:contrast +pipe:in
  !cov [contrast_formula_spec]

@s3m contrast_weights.contrast_set (x:s3:contrast_set, term:s3:event_term, ...) | Compute contrast weights for set -> lst
  +family:contrast +pipe:in
  !cov [contrast_set]

@s3m contrast_weights.event_model (x:s3:event_model, ...) | Compute contrast weights for event model -> lst
  +family:contrast +pipe:in
  !cov [event_model]

@s3m contrast_weights.oneway_contrast_spec (x:s3:oneway_contrast_spec, term:s3:event_term, ...) | Compute one-way contrast weights -> lst
  +family:contrast +pipe:in
  !cov [oneway_contrast_spec]

@s3m contrast_weights.pair_contrast_spec (x:s3:pair_contrast_spec, term:s3:event_term, ...) | Compute pair contrast weights -> lst
  +family:contrast +pipe:in
  !cov [pair_contrast_spec]

@s3m contrast_weights.poly_contrast_spec (x:s3:poly_contrast_spec, term:s3:event_term, ...) | Compute polynomial contrast weights -> lst
  +family:contrast +pipe:in
  !cov [poly_contrast_spec]

@s3m contrast_weights.unit_contrast_spec (x:s3:unit_contrast_spec, term:s3:event_term, ...) | Compute unit contrast weights -> lst
  +family:contrast +pipe:in
  !cov [unit_contrast_spec]

@s3m design_matrix.baseline_model (x:s3:baseline_model, blockid:int?=NULL, allrows:lgl?=FALSE, ...) | Extract design matrix -> df
  +family:design +pipe:in
  !cov [baseline_model]

@s3m design_matrix.event_model (x:s3:event_model, blockid:int?=NULL, ...) | Extract design matrix -> df
  +family:design +pipe:in
  !cov [event_model]

@s3m design_matrix.event_term (x:s3:event_term, drop.empty:lgl?=TRUE, ...) | Extract design matrix -> df
  +family:design +pipe:in
  !cov [event_term]

@s3m plot.baseline_model (x:s3:baseline_model, term_name:chr?=NULL, title:chr?=NULL, xlab:chr?="Time", ylab:chr?="Design Matrix Value", line_size:dbl?=1, color_palette:chr?="Set1", ...) | Plot baseline model -> Gg
  +family:plot +pipe:in +side:plot[device]
  !cov [baseline_model]

@s3m plot.event_model (x:s3:event_model, term_name:chr?=NULL, ...) | Plot event model -> Gg
  +family:plot +pipe:in +side:plot[device]
  !cov [event_model]

@s3m plot_contrasts.event_model (x:s3:event_model, absolute_limits:lgl?=FALSE, rotate_x_text:lgl?=TRUE, scale_mode:chr["auto"|"diverging"|"one_sided"]?="auto", coord_fixed:lgl?=TRUE, ...) | Plot contrasts -> Gg
  +family:plot +pipe:in +side:plot[device]
  !cov [event_model]

@s3m print.baseline_model (x:s3:baseline_model, ...) | Print baseline model -> invisible
  +family:print +pipe:in +side:console[print]
  !cov [baseline_model]

@s3m print.contrast (x:s3:contrast, ...) | Print contrast -> invisible
  +family:print +pipe:in +side:console[print]
  !cov [contrast]

@s3m print.contrast_set (x:s3:contrast_set, ...) | Print contrast set -> invisible
  +family:print +pipe:in +side:console[print]
  !cov [contrast_set]

@s3m print.event_model (x:s3:event_model, ...) | Print event model -> invisible
  +family:print +pipe:in +side:console[print]
  !cov [event_model]

@s3m print.event_term (x:s3:event_term, ...) | Print event term -> invisible
  +family:print +pipe:in +side:console[print]
  !cov [event_term]

@s3m regressors.event_term (x:s3:event_term, hrf:s3:HRF, sampling_frame:s3:sampling_frame, summate:lgl?=FALSE, drop.empty:lgl?=TRUE, ...) | Extract regressors -> lst
  +family:regressor +pipe:in
  !cov [event_term]

@s3m split_onsets.event_term (x:s3:event_term, sframe:s3:sampling_frame, global:lgl?=FALSE, blocksplit:lgl?=FALSE, ...) | Split onsets -> lst<vec<dbl>>
  +family:onset +pipe:in
  !cov [event_term]

## Dependencies
- Imports: fmrihrf, stats, assertthat, rlang, stringr, dplyr, tidyr, purrr, tibble, Matrix, splines, plotly, ggplot2, utils, cli
- Suggests: testthat, knitr, rmarkdown, covr

---
## Meta-Footer
- Micro-DSL Version: v3.4-source
- Package: fmridesign (Version: 0.1.0)
- Generated: 2023-10-07T12:00:00Z
- Features: types[constrained] sigils[specific] metadata[rich] semantics[annotated]
- Coverage: 52 / 52 exports
- Provenance: exports[NAMESPACE], enums[match.arg/switch], constraints[assertions/checks]
```