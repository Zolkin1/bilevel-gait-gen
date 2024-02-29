//
// Created by zolkin on 2/28/24.
//

#include "optitrack-stream-client/client.h"

int main() {
    stream_client::Client client(300);
    std::cout << "[Optitrack Client] initializing listening." << std::endl;
    client.Listen();
}