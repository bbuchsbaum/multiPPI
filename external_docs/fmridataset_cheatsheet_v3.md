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
@pkg fmridataset | Unified Container for fMRI Datasets

## Type Aliases:
DF = df|tbl|data.table
V<T> = vec<T>
Fml = s3:formula
Gg = s3:ggplot

## Constraint Definitions:
@constraint positive_weights | Positive numeric weights
  type: vec<dbl>
  validates: all(. > 0)

# fmridataset

## 1. Core

@f fmri_dataset (scans:chr[], mask:chr, TR:dbl, run_length:int[], event_table:DF?=data.frame(), base_path:chr?=".", censor:vec<int>?=NULL, preload:lgl?=FALSE, mode:chr["normal"|"bigvec"|"mmap"|"filebacked"]?="normal", dummy_mode:lgl?=FALSE) | Create fMRI dataset -> s3:fmri_dataset
  +family:dataset +pipe:in +side:fs[read] +see:fmri_mem_dataset, fmri_h5_dataset
  - scans : @requires:files-exist @purpose:data-files
  - mask : @requires:file-exist @purpose:mask-file
  - TR : @units:seconds @purpose:temporal-resolution
  - run_length : @purpose:run-lengths
  - event_table : @purpose:event-metadata
  - base_path : @purpose:relative-paths
  - censor : @purpose:exclude-scans
  - preload : @purpose:load-strategy
  - mode : @purpose:storage-mode
  - dummy_mode : @purpose:test-mode

@f fmri_mem_dataset (scans:lst<s3:NeuroVec>, mask:s3:NeuroVol, TR:dbl, run_length:int[], event_table:DF?=data.frame(), base_path:chr?=".", censor:vec<int>?=NULL) | Create in-memory fMRI dataset -> s3:fmri_mem_dataset
  +family:dataset +pipe:in +side:fs[read] +see:fmri_dataset, fmri_h5_dataset
  - scans : @requires:NeuroVec-objects @purpose:data-objects
  - mask : @requires:NeuroVol-object @purpose:mask-object
  - TR : @units:seconds @purpose:temporal-resolution
  - run_length : @purpose:run-lengths
  - event_table : @purpose:event-metadata
  - base_path : @purpose:relative-paths
  - censor : @purpose:exclude-scans

@f fmri_h5_dataset (h5_files:chr[], mask_source:chr|s3:NeuroVol, TR:dbl, run_length:int[], event_table:DF?=data.frame(), base_path:chr?=".", censor:vec<int>?=NULL, preload:lgl?=FALSE, mask_dataset:chr?="data/elements", data_dataset:chr?="data") | Create fMRI dataset from H5 files -> s3:fmri_file_dataset
  +family:dataset +pipe:in +side:fs[read] +see:fmri_dataset, fmri_mem_dataset
  - h5_files : @requires:files-exist @purpose:data-files
  - mask_source : @requires:file-exist @purpose:mask-file
  - TR : @units:seconds @purpose:temporal-resolution
  - run_length : @purpose:run-lengths
  - event_table : @purpose:event-metadata
  - base_path : @purpose:relative-paths
  - censor : @purpose:exclude-scans
  - preload : @purpose:load-strategy
  - mask_dataset : @purpose:mask-path
  - data_dataset : @purpose:data-path

@f latent_dataset (source:chr[]|lst<s3:LatentNeuroVec>, TR:dbl, run_length:int[], event_table:DF?=data.frame(), base_path:chr?=".", censor:vec<int>?=NULL, preload:lgl?=FALSE) | Create latent space fMRI dataset -> s3:latent_dataset
  +family:dataset +pipe:in +side:fs[read] +see:fmri_dataset, fmri_mem_dataset
  - source : @requires:files-exist @purpose:data-files
  - TR : @units:seconds @purpose:temporal-resolution
  - run_length : @purpose:run-lengths
  - event_table : @purpose:event-metadata
  - base_path : @purpose:relative-paths
  - censor : @purpose:exclude-scans
  - preload : @purpose:load-strategy

@f get_latent_scores (x:s3:latent_dataset, rows:int[]?=NULL, cols:int[]?=NULL, ...:expr) | Extract latent scores -> mat<dbl>
  +family:data-access +pipe:in +see:get_spatial_loadings, reconstruct_voxels
  - x : @purpose:dataset
  - rows : @purpose:subset-rows
  - cols : @purpose:subset-cols

@f get_spatial_loadings (x:s3:latent_dataset, components:int[]?=NULL, ...:expr) | Extract spatial loadings -> mat<dbl>
  +family:data-access +pipe:in +see:get_latent_scores, reconstruct_voxels
  - x : @purpose:dataset
  - components : @purpose:subset-components

@f get_component_info (x:s3:latent_dataset, ...:expr) | Get component metadata -> lst
  +family:data-access +pipe:in +see:get_latent_scores, get_spatial_loadings
  - x : @purpose:dataset

@f reconstruct_voxels (x:s3:latent_dataset, rows:int[]?=NULL, voxels:int[]?=NULL, ...:expr) | Reconstruct voxel data -> mat<dbl>
  +family:data-access +pipe:in +see:get_latent_scores, get_spatial_loadings
  - x : @purpose:dataset
  - rows : @purpose:subset-rows
  - voxels : @purpose:subset-voxels

@f fmri_series (dataset:s3:fmri_dataset, selector:any?=NULL, timepoints:int[]?=NULL, output:chr["fmri_series"|"DelayedMatrix"]?="fmri_series", event_window:any?=NULL, ...:expr) | Query fMRI time series -> s3:fmri_series|DelayedMatrix
  +family:data-access +pipe:in +see:fmri_dataset, fmri_mem_dataset
  - dataset : @purpose:dataset
  - selector : @purpose:spatial-selector
  - timepoints : @purpose:temporal-selector
  - output : @purpose:return-format

@f fmri_group (subjects:DF, id:chr, dataset_col:chr?="dataset", space:chr?=NULL, mask_strategy:chr["subject_specific"|"intersect"|"union"]?="subject_specific") | Create fMRI group -> s3:fmri_group
  +family:group +pipe:in +see:fmri_dataset, fmri_mem_dataset
  - subjects : @purpose:subject-data
  - id : @purpose:identifier
  - dataset_col : @purpose:dataset-column
  - space : @purpose:common-space
  - mask_strategy : @purpose:mask-strategy

@f iter_subjects (gd:s3:fmri_group, order_by:chr?=NULL) | Iterate subjects -> lst
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - order_by : @purpose:order-iteration

@f group_map (gd:s3:fmri_group, .f:fml, ..., out:chr["list"|"bind_rows"]?="list", order_by:chr?=NULL, on_error:chr["stop"|"warn"|"skip"]?="stop") | Map function over subjects -> lst|DF
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - .f : @purpose:map-function
  - ... : @purpose:additional-args
  - out : @purpose:output-format
  - order_by : @purpose:order-iteration
  - on_error : @purpose:error-handling

@f group_reduce (gd:s3:fmri_group, .map:fml, .reduce:fml, .init:any, order_by:chr?=NULL, on_error:chr["stop"|"warn"|"skip"]?="stop", ...) | Reduce over subjects -> any
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - .map : @purpose:map-function
  - .reduce : @purpose:reduce-function
  - .init : @purpose:initial-value
  - order_by : @purpose:order-iteration
  - on_error : @purpose:error-handling

@f filter_subjects (gd:s3:fmri_group, ...:expr) | Filter subjects -> s3:fmri_group
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - ... : @purpose:filter-expressions

@f mutate_subjects (gd:s3:fmri_group, ...:expr) | Mutate subject attributes -> s3:fmri_group
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - ... : @purpose:mutate-expressions

@f left_join_subjects (gd:s3:fmri_group, y:DF, by:chr?=NULL, ...) | Left join subject metadata -> s3:fmri_group
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - y : @purpose:join-data
  - by : @purpose:join-keys

@f sample_subjects (gd:s3:fmri_group, n:int, replace:lgl?=FALSE, strata:chr?=NULL) | Sample subjects -> s3:fmri_group
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - n : @purpose:sample-size
  - replace : @purpose:sample-replacement
  - strata : @purpose:sample-strata

@f stream_subjects (gd:s3:fmri_group, prefetch:int?=1, order_by:chr?=NULL) | Stream subjects -> lst
  +family:group +pipe:in +see:fmri_group, fmri_dataset
  - gd : @purpose:group-dataset
  - prefetch : @purpose:prefetch-size
  - order_by : @purpose:order-iteration

## Dependencies
- Imports: assertthat, cachem, deflist, fmrihrf, fs, lifecycle, memoise, Matrix, methods, neuroim2, purrr, tibble, DelayedArray, S4Vectors, utils
- Suggests: bench, bidser, crayon, arrow, dplyr, fmristore, foreach, mockery, mockr, Rarr, S4Arrays, testthat (>= 3.0.0), knitr, rmarkdown

---
## Meta-Footer
- Micro-DSL Version: v3.4-source
- Package: fmridataset (Version: 0.8.9)
- Generated: 2023-10-05T12:00:00Z
- Features: types[constrained] sigils[specific] metadata[rich] semantics[annotated]
- Coverage: 50 / 50 exports
- Provenance: exports[NAMESPACE], enums[match.arg/switch], constraints[assertions/checks]
```