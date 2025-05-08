#include "pregel_app_mst.h"

int main(int argc, char *argv[])
{
    init_workers();
    pregel_mst("/toyFolder", "/toyOutput");
    worker_finalize();
    return 0;
}
