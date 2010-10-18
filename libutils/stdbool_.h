
#ifndef _LIBUTILS_STDBOOL_H
#define _LIBUTILS_STDBOOL_H

typedef enum { must_promote_to_int = -1, false = 0, true = 1 } utils_bool;

#define bool utils_bool
#define false 0
#define true 1
#define __bool_true_false_are_defined 1

#endif
