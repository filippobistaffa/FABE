#ifndef MAIN_H_
#define MAIN_H_

#include "types.hpp"
#include <limits>       // std::numeric_limits

// known threshold values to remove rows
static const char* datasets[] = { "spot5",
                                  "celar6",
                                  "mastermind",
                                  "iscas89",
                                  "uai"
                                };

static const value thresholds[] = { 100000,
                                    1000,
                                    100,
                                    100,
                                    std::numeric_limits<value>::infinity()
                                  };

#define N_DATASETS (sizeof(datasets) / sizeof(char*))

#endif /* MAIN_H_ */
