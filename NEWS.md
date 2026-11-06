# flanders 2.3.11
## Bug fixes
* `anndata2coloc_row()` no longer crashes with "invalid regular expression"
  when `t1_phenotype_id` / `t2_phenotype_id` contains characters such as
  `:`, `;`, `,` or is unusually long (e.g. an aggregated list of cell
  types and regions).
## Output schema
* `anndata2coloc()` / `anndata2coloc_row()` Input and output colocalization
  tables carry the same columns.
