!*****************************************************************************************
!
!  The main module that uses all the other modules.
!  Allows for a single `use seiscope_optimization_toolbox`
!  to access the entire library.

module seiscope_optimization_toolbox

use typedef
use miscellaneous
use opt_PSTD
use opt_PNLCG
use opt_LBFGS
use opt_PLBFGS
use opt_TRN
use opt_PTRN

implicit none

public

end module seiscope_optimization_toolbox
!*****************************************************************************************
