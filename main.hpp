#ifndef MAIN_H_
#define MAIN_H_

#include "io.hpp"
#include "order.hpp"
#include "conversion.hpp"
#include "be.hpp"

#define PROFILE "trace.prof"

#ifdef PROFILE
#include <gperftools/profiler.h>
#endif

// known threshold values to remove rows
static const char* datasets[] = { "spot5",
                                  "mastermind",
                                  "iscas89"
                                };

static const value thresholds[] = { 100000,
                                    100,
                                    100
                                  };

#define N_DATASETS (sizeof(datasets) / sizeof(char*))

#endif /* MAIN_H_ */
