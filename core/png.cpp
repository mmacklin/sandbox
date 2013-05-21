#include "png.h"
#include "log.h"

#include <iostream>

#include "external/stb_image/stb_image.c"

bool PngLoad(const char* filename, PngImage& image)
{	
	int x, y, c;
	
	uint8_t* data = stbi_png_load(filename, &x, &y, &c, 4);
	
	if (data)
	{
		int s = x*y;
		
		image.m_data = new uint32_t[s];
		memcpy(image.m_data, data, s*4);
		
		image.m_width = (uint16)x;
		image.m_height = (uint16)y;
		
		stbi_image_free(data);
		
		Log::Info << "Loaded png " << filename << image.m_width << ", " << image.m_height << std::endl;		
		return true;
	}
	else 
	{
		Log::Info << "Could not find " << filename << " for loading" << std::endl;
		return false;
	}
}

#if 0

bool PngLoad(const char* filename, PngImage& image)
{
	png_structp png_ptr;
	png_infop info_ptr;
	//unsigned int sig_read = 0;
	png_uint_32 width, height;
	int bit_depth, color_type, interlace_type;
	FILE *fp;

	if ((fp = fopen(filename, "rb")) == NULL)
		return NULL;

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
		NULL, NULL, NULL);

	if (png_ptr == NULL)
	{
		fclose(fp);
		return NULL;
	}

	/* Allocate/initialize the memory for image information.  REQUIRED. */
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL)
	{
		fclose(fp);
		png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
		return false;
	}

	/* Set error handling if you are using the setjmp/longjmp method (this is
	* the normal method of doing things with libpng).  REQUIRED unless you
	* set up your own error handlers in the png_create_read_struct() earlier.
	*/
	if (setjmp(png_ptr->jmpbuf))
	{
		/* Free all of the memory associated with the png_ptr and info_ptr */
		png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
		fclose(fp);
		/* If we get here, we had a problem reading the file */
		return false;
	}

	png_init_io(png_ptr, fp);
	png_read_info(png_ptr, info_ptr);
	png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
		&interlace_type, NULL, NULL);
	
	/* Add filler (or alpha) uint8_t (before/after each RGB triplet) */
	png_set_expand(png_ptr);
	png_set_filler(png_ptr, 0xff, PNG_FILLER_AFTER);
	//png_set_gray_1_2_4_to_8(png_ptr);
	png_set_palette_to_rgb(png_ptr);
	png_set_gray_to_rgb(png_ptr);
	
	//png_set_bgr(png_ptr);

	int aNumBytes = png_get_rowbytes(png_ptr, info_ptr) * height;

	uint32_t* aBits = new uint32_t[aNumBytes/4];

	for (uint32_t i = 0; i < height; i++)
	{
		png_bytep anAddr = (png_bytep) &aBits[i*width];
		png_read_rows(png_ptr, (png_bytepp) &anAddr, NULL, 1);
	}
	
	std::cout << filename << " " << width << " " << height << " : " << interlace_type << std::endl;
	
	image.m_data = (uint32_t*)aBits;
	image.m_width = (uint16)width;
	image.m_height = (uint16)height;

	return true;
}

#endif
