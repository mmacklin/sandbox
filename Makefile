CTAGS  = ctags --fields=+l --c-kinds=+p --c++-kinds=+p 

FILES  = core/*.h 
FILES += external/glut/glut.h 
FILES += /Library/Frameworks/CUDA.framework/Versions/A/Headers/*.h
FILES += /usr/local/cuda/include/cuda_runtime_api.h
FILES += /Developer/SDKs/MacOSX10.6.sdk/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers/gl.h 

tags: $(FILES) Makefile
	$(CTAGS) $(FILES) 
	
