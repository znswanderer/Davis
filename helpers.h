#ifndef HELPERS_H
#define HELPERS_H

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// Attention: The mod (%) operator in C works 
// differently in the negative domain from Python:
// python: -1 % 4 = 3
// C:      -1 % 4 = -1
// We use the TRUNCED macro to get python's behaviour.
// Please note: this version only works for x >= -max_x
#define TRUNCED(x, max_x) (((x) + (max_x)) % (max_x))


#endif
