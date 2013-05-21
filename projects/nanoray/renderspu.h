#ifdef PLATFORM_PS3

#include <libspe2.h>

extern spe_program_handle_t gRenderJobSpuModule;

// on SPU we need to kick create an SPU job context and wait for that
void* RenderJobThreadFuncSPU(void* param)
{
	int retval;
	unsigned int entry_point = SPE_DEFAULT_ENTRY; 
	spe_context_ptr_t my_context;

	// create the SPE Context 
	my_context = spe_context_create(0, 0);

	// load the embedded code into this context
	spe_program_load(my_context, &gRenderJobSpuModule);

	// Run the SPE program until completion
	do {
		retval = spe_context_run(my_context, &entry_point, 0, param, NULL, NULL);
	} while (retval > 0); // Run until exit or error 

	pthread_exit(NULL);
}

#endif