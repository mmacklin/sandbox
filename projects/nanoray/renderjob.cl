/*!
 * Sample kernel which multiplies every element of the input array with
 * a constant and stores it at the corresponding output array
 */

 #define PLATFORM_OPENCL 1

 #include "renderjob.h"


__kernel void renderKernel(__global float4* output)
{
    uint tid = get_global_id(0);
    
    output[tid] =  (float4)(tid*0.00001f, tid*0.00001f, 0.0f, 1.0f);
}
