# Generalizability Theory ICC toolbox

*Purpose:* Calculate ICC coefficients including multiple sources of error (i.e., facets) and perform Decision Study.

## Getting Started

### Prerequisites

Matlab (tested with Matlab 2018b).

### Usage

1. Create data and factor table variables. This can be accomplished using the *load\_reliability\_data* script, which automatically creates a factor table from the filenames. This can be changed for your data or data can be loaded manually (try "help load\_reliability\_data" for more info on the format). 

2. Run reliability analysis. Example: *[icc,var,stats,sigmask]=run\_reliability(data',data,'factors',ftbl)*. Try "help run\_reliability" for more info. This will result in average D- and G-coefficients, a Decision Study visualization, and a variable "icc" which contains detailed ICC results.

*References:* 

Noble, S., Spann, M. N., Tokoglu, F., Shen, X., Constable, R. T., & Scheinost, D. (2017). Influences on the test–retest reliability of functional connectivity MRI and its relationship with behavioral utility. Cerebral Cortex, 27(11), 5415-5429.

Noble, S., Scheinost, D., Finn, E. S., Shen, X., Papademetris, X., McEwen, S. C., ... & Mirzakhanian, H. (2017). Multisite reliability of MR-based functional connectivity. Neuroimage, 146, 959-970.

Webb, N.M., Shavelson, R.J. and Haertel, E.H., 2006. Reliability coefficients and generalizability theory. Handbook of statistics, 26, pp.81-124.
