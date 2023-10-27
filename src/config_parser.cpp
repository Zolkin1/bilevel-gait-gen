//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"

namespace mpc::utils {
    ConfigParser::ConfigParser(const std::string& file_name) : file_name_(file_name) {
        config_ = YAML::LoadFile(file_name_);
    }

    Eigen::VectorXd ConfigParser::ParseEigenVector(const std::string& element) {
        auto temp = config_[element].as<std::vector<double>>();
        Eigen::VectorXd vec(temp.size());
        for (int i = 0; i < temp.size(); i++) {
            vec(i) = temp.at(i);
        }

        return vec;
    }

    double ConfigParser::ParseNumber(const std::string &element) {
        return config_[element].as<double>();
    }

    std::string ConfigParser::ParseString(const std::string &element) {
        return config_[element].as<std::string>();
    }

    std::vector<std::string> ConfigParser::ParseStringVector(const std::string &element) {
        return config_[element].as<std::vector<std::string>>();
    }

    std::vector<int> ConfigParser::ParseIntVector(const std::string &element) {
        return config_[element].as<std::vector<int>>();
    }
} //mpc::utils