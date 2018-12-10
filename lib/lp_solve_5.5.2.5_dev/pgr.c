#include <stdio.h>
#include "test/lp_lib.h"
main()
{
lprec *lp;
lp=make_lp(0,4);
delete_lp(lp);
}
