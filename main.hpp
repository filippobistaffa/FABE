#ifndef MAIN_H_
#define MAIN_H_

#include "io.hpp"
#include "order.hpp"
#include "conversion.hpp"
#include "be.hpp"

#include <chrono>
#include <numeric>      // accumulate

// known threshold values to remove rows
static const char* datasets[] = { "spot5",
                                  "celar6",
                                  "mastermind",
                                  "iscas89"
                                };

static const value thresholds[] = { 100000,
                                    1000,
                                    100,
                                    100
                                  };

#define N_DATASETS (sizeof(datasets) / sizeof(char*))

#endif /* MAIN_H_ */
