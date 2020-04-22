#ifndef LOG_HPP_
#define LOG_HPP_

#include <iostream>     // cout
#include <iomanip>      // setw

using namespace std;

#define TOTAL_WIDTH 79
#define COLUMN_WIDTH ((TOTAL_WIDTH - 7) / 2)

static float progress;

__attribute__((always_inline)) inline
void log_line() {

        cout << "+";
        for (size_t i = 0; i < COLUMN_WIDTH + 2; ++i) {
                cout << "-";
        }
        cout << "+";
        for (size_t i = 0; i < COLUMN_WIDTH + 2; ++i) {
                cout << "-";
        }
        cout << "+" << endl;
}

template <typename T>
__attribute__((always_inline)) inline
void log_value(string name, T val, bool filename = false) {

        cout << "| ";
        if (name.length() > COLUMN_WIDTH) {
                cout << name.substr(0, COLUMN_WIDTH - 3) << "...";
        } else {
                cout << setw(COLUMN_WIDTH) << left << name;
        }
        cout << " | ";
        if constexpr (is_same<T, char *>::value) {
                const string str = string(val);
                if (filename) {
                        cout << "..." << str.substr(str.length() - COLUMN_WIDTH + 3);
                } else {
                        cout << str.substr(0, COLUMN_WIDTH - 3) << "...";
                }
        } else {
                cout << setw(COLUMN_WIDTH) << left << val;
        }
        cout << " |" << endl;
}

__attribute__((always_inline)) inline
void log_progress_increase(float step, float tot) {

        if (progress == tot) {
                return;
        }

        if (progress == 0) {
                cout << "|";
        }

        const size_t cur_p = (TOTAL_WIDTH - 2) * (progress / tot);
        const size_t new_p = (TOTAL_WIDTH - 2) * ((progress + step) / tot);

        for (size_t i = cur_p; i < new_p; ++i) {
                cout << "·" << flush;;
        }

        progress += step;

        if (progress == tot) {
                cout << "|" << endl;
        }
}

#endif /* LOG_HPP_ */