//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CONFIG_PARSER_H
#define BILEVEL_GAIT_GEN_CONFIG_PARSER_H

#include <Eigen/Core>

#include "yaml-cpp/yaml.h"

namespace utils {
    class ConfigParser {
    public:
        ConfigParser(const std::string& file_name);

        Eigen::VectorXd ParseEigenVector(const std::string& element);

        template<typename scalar>
        scalar ParseNumber(const std::string& element) {
            return config_[element].as<scalar>();
        }

        std::string ParseString(const std::string& element);

        // TODO: Remove
        std::vector<std::string> ParseStringVector(const std::string& element);

        template<typename scalar>
        std::vector<scalar> ParseStdVector(const std::string& element) {
            return config_[element].as<std::vector<scalar>>();
        }

    private:
        std::string file_name_;

        YAML::Node config_;
    };
} // mpc::utils


#endif //BILEVEL_GAIT_GEN_CONFIG_PARSER_H
