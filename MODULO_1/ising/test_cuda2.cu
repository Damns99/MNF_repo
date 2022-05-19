#include <cstdio>

class myClass
{
public:
        bool bool_var;    // Set from host and readable from device
        int  data_size;   // Set from host
        __host__ myClass();
        __host__ ~myClass();
        __host__ void setValues(bool iftrue, int size);
        __device__ void dosomething(int device_parameter);
        __host__ void export_data();

        // completely unknown methods
        __host__ void prepareDeviceObj();
        __host__ void retrieveDataToHost();
        __host__ void runKernel();
private:
        int *data; // Filled in device, shared between threads, at the end copied back to host for data output
        int *h_data;
};

__host__ myClass::myClass()
{
}

__host__ myClass::~myClass()
{
}

__host__ void myClass::prepareDeviceObj(){
        cudaMemcpy(data, h_data, data_size*sizeof(h_data[0]), cudaMemcpyHostToDevice);
}
__host__ void myClass::retrieveDataToHost(){
        cudaMemcpy(h_data, data, data_size*sizeof(h_data[0]), cudaMemcpyDeviceToHost);
}

__host__ void myClass::setValues(bool iftrue, int size)
{
        bool_var  = iftrue;
        data_size = size;
        cudaMalloc(&data, data_size*sizeof(data[0]));
        h_data = (int *)malloc(data_size*sizeof(h_data[0]));
        memset(h_data, 0, data_size*sizeof(h_data[0]));
}

__device__ void myClass::dosomething(int idx)
{
        int toadd = idx+data_size;
        atomicAdd(&(data[idx]), toadd); // data should be unique among threads
}
__host__ void myClass::export_data(){
        for (int i = 0; i < data_size; i++) printf("%d ", h_data[i]);
        printf("\n");
        cudaFree(data);
        free(h_data);
}


__global__ void myKernel(myClass obj)
{
        const int idx = blockIdx.x*blockDim.x + threadIdx.x;
        if(idx < obj.data_size)
        {
                if(!obj.bool_var)
                        printf("Object is not up to any task here!");
                else
                {
                        //printf("Object is ready!");
                        obj.dosomething(idx);
                }
        }
}

__host__ void myClass::runKernel() {
    int some_number = 40;
    myKernel<<<1,some_number>>>(*this);
}


myClass globalInstance;

int main(int argc, char** argv)
{
        int some_number = 40;
        globalInstance.setValues(true, some_number);
        globalInstance.prepareDeviceObj();
        globalInstance.runKernel();
        globalInstance.retrieveDataToHost();
        globalInstance.export_data();
        exit(EXIT_SUCCESS);
}