#include "rendercl.h"
#include "renderjob.h"
#include "renderjobqueue.h"

#include "shared/platform.h"

#include <cl/cl.h>
#include <iostream>
#include <string>

#if _DEBUG
#define CL_CHECK_ERROR(x) { assert(x == CL_SUCCESS); }
#else
#define CL_CHECK_ERROR(x) x
#endif

inline void build(cl_program program, cl_device_id* devices, const char* options = "")
{
	// build
	cl_int ret_val = clBuildProgram(program, 1, devices, options, NULL, NULL);

	// avoid abortion due to CL_BILD_PROGRAM_FAILURE
	if (ret_val != CL_SUCCESS && ret_val != CL_BUILD_PROGRAM_FAILURE)
		CL_CHECK_ERROR(ret_val);

	cl_build_status build_status;
	CL_CHECK_ERROR(clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL));
	if (build_status == CL_SUCCESS)
 		return;

	char *build_log;
	size_t ret_val_size;
	CL_CHECK_ERROR(clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size));

	build_log = new char[ret_val_size+1];
	CL_CHECK_ERROR(clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL));

	// to be carefully, terminate with \0
	// there's no information in the reference whether the string is 0 terminated or not
	build_log[ret_val_size] = '\0';

	std::cout << build_log << std::endl;

	delete[] build_log;
}

using namespace std;

void* RenderJobThreadFuncOpenCL(void* param)
{
	RenderJobQueue* queue = (RenderJobQueue*)(param);

	// pop off render tiles until there are none left
	RenderJob* job = queue->Pop();
	if (!job)
		return 0;

	cl_int status = 0;

    cl_device_type dType = CL_DEVICE_TYPE_GPU;
    cl_context context = clCreateContextFromType(0, dType, NULL, NULL, &status);

	if (status != CL_SUCCESS)
	{
		cout << "clCreateContextFromType failed." << endl;
        return 0;
	}
	
	size_t contextDescriptorSize;
	status = clGetContextInfo(context, CL_CONTEXT_DEVICES,0, 0, &contextDescriptorSize);

	if (status != CL_SUCCESS)
	{
		cout << "clCreateContextFromType failed." << endl;
        return 0;
	}

	cl_device_id* devices = (cl_device_id*)malloc(contextDescriptorSize);
	status = clGetContextInfo(context, CL_CONTEXT_DEVICES, contextDescriptorSize, devices, 0);

	if (status != CL_SUCCESS)
	{
		cout << "clGetContextInfo failed." << endl;
        return 0;
	}
	
	// some stats
	size_t paramAddressBits;
	clGetDeviceInfo(devices[0], CL_DEVICE_ADDRESS_BITS, sizeof(paramAddressBits), &paramAddressBits, NULL);
	cout << "CL_DEVICE_ADDRESS_BITS: " << paramAddressBits << endl;

	cl_command_queue cmdQueue;
	cmdQueue = clCreateCommandQueue(context, devices[0], 0, 0);
	
	// create cl buffer for render job output
	const uint32_t bufSize = job->m_renderRect.Width() * job->m_renderRect.Height() * 4 * sizeof(cl_float);

    cl_mem outputBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, bufSize, NULL, &status);

	if (status != CL_SUCCESS)
	{
		cout << "clCreateBuffer failed." << endl;
        return 0;
	}
	
	string sfile = LoadFileToString("RenderJob.cl");

	const char* source = sfile.c_str();
	size_t sourceSize[] = { sfile.length() };
	
	cl_program program = clCreateProgramWithSource(context, 1, &source, sourceSize, &status);

	if (status != CL_SUCCESS)
	{
		cout << "clCreateProgramWithSource failed." << endl;
        return 0;
	}

	// build the program
	build(program, devices, "-I.");

    // get a kernel object handle for a kernel with the given name 
    cl_kernel kernel = clCreateKernel(program, "renderKernel", &status);

	if (status != CL_SUCCESS)
	{
		cout << "clCreateKernel failed." << endl;
        return 0;
	}

	// set kernel params
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&outputBuffer);
	
	if (status != CL_SUCCESS)
	{
		cout << "clSetKernelArg failed." << endl;
        return 0;
	}

	size_t globalThreads[1];
    size_t localThreads[1];
	cl_event events[2];

	globalThreads[0] = job->m_renderRect.Width()*job->m_renderRect.Height();
    localThreads[0]  = 1;

	status = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, globalThreads, localThreads, 0, NULL, &events[0]);

	if(status != CL_SUCCESS) 
	{ 
		cout << "clEnqueueNDRangeKernel failed" << endl;
		return 0;
	}

    status = clWaitForEvents(1, &events[0]);
	if(status != CL_SUCCESS) 
	{ 
		cout << "clWaitForEvents failed" << endl;
		return 0;
	}	

	clReleaseEvent(events[0]);

	// read buffer back
	status = clEnqueueReadBuffer(cmdQueue, outputBuffer, CL_TRUE, 0, bufSize, job->m_output, 0, 0, 0);
	
	if(status != CL_SUCCESS) 
	{ 
		cout << "clEnqueueReadBuffer failed" << endl;
		return 0;
	}	

	Sleep(2.0);

    return NULL;
}

