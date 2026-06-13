/**************************************************************************
 * pylogging.hpp
 *
 * Copyright 2026 Gabriel Machado
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **************************************************************************/

#pragma once

#include <pybind11/pybind11.h>


class PyLogger {
private:
    py::object _logger;
    py::object _log_handle;
    py::object _log_make_record;
    bool _is_enabled_cache[5] = {false};
    static inline PyLogger *_global_logger = nullptr;

    inline void _set_logger(py::object logger) {
        py::gil_scoped_acquire acquire;
        _logger = logger;
        _log_make_record = _logger.attr("makeRecord");
        _log_handle = _logger.attr("handle");

        auto is_enabled_for = _logger.attr("isEnabledFor");
        for (int i = 0; i < 5; ++i) {
            auto is_enabled = is_enabled_for(i * 10);
            _is_enabled_cache[i] = is_enabled.cast<bool>();
        }
    }

public:
    PyLogger() {
        py::gil_scoped_acquire acquire;
        auto logging = py::module::import("logging");
        _set_logger(logging.attr("getLogger")());

        py_assert(!_global_logger, "A global logger already exists");
        _global_logger = this;
    }

    ~PyLogger() {
        py::gil_scoped_acquire acquire;
        _global_logger = nullptr;
    }

    static inline void set_logger(py::object logger) {
        if (logger.is_none()) {
            delete _global_logger;
            _global_logger = nullptr;
            return;
        }

        if (_global_logger == nullptr)
            new PyLogger();

        _global_logger->_set_logger(logger);
    }

    inline py::object get_logger() const {
        py::gil_scoped_acquire acquire;
        return _logger;
    }

    inline void pylog(int level, const char *file, int lineno, const char *msg, py::tuple&& args, const char *func) {
        if (!_is_enabled_cache[level / 10])
            return;

        _log_handle(_log_make_record("par.rs", level, file, lineno, msg, args, py::none(), func));
    }

    static inline PyLogger *instance() {
        return _global_logger;
    }
};


#define _pylog(level, msg, ...) if (PyLogger::instance()) PyLogger::instance()->pylog(level, __FILE__, __LINE__, msg, py::make_tuple(__VA_ARGS__), __FUNCTION__)
#define log_critical(msg, ...) _pylog(50, msg, __VA_ARGS__)
#define log_error(msg, ...) _pylog(40, msg, __VA_ARGS__)
#define log_warning(msg, ...) _pylog(30, msg, __VA_ARGS__)
#define log_info(msg, ...) _pylog(20, msg, __VA_ARGS__)
#define log_debug(msg, ...) _pylog(10, msg, __VA_ARGS__)
