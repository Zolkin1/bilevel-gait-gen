//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CONFIG_PARSER_H
#define BILEVEL_GAIT_GEN_CONFIG_PARSER_H

#include <Eigen/Core>

#include "yaml-cpp/yaml.h"

namespace mpc::utils {
    class ConfigParser {
    public:
        ConfigParser(const std::string& file_name);

        Eigen::VectorXd ParseEigenVector(const std::string& element);

        double ParseNumber(const std::string& element);

        std::string ParseString(const std::string& element);

    private:
        std::string file_name_;

        YAML::Node config_;
    };
} // mpc::utils


#endif //BILEVEL_GAIT_GEN_CONFIG_PARSER_H
