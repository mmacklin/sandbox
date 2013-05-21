#include "maths.h"

uint32_t seed1;
uint32_t seed2;

void RandInit()
{
	seed1 = 315645664;
	seed2 = seed1 ^ 0x13ab45fe;
}

static float s_identity[4][4] = { {1.0f, 0.0f, 0.0f, 0.0f},
								  {0.0f, 1.0f, 0.0f, 0.0f},
								  {0.0f, 0.0f, 1.0f, 0.0f},
								  {0.0f, 0.0f, 0.0f, 1.0f} };

template <>
XMatrix44<float> XMatrix44<float>::kIdentity(s_identity[0]);
								

Colour::Colour(Colour::Preset p)
{
	switch (p)
	{
	case kRed:
		*this = Colour(1.0f, 0.0f, 0.0f);
		break;
	case kGreen:
		*this = Colour(0.0f, 1.0f, 0.0f);
		break;
	case kBlue:
		*this = Colour(0.0f, 0.0f, 1.0f);
		break;
	case kWhite:
		*this = Colour(1.0f, 1.0f, 1.0f);
		break;
	case kBlack:
		*this = Colour(0.0f, 0.0f, 0.0f);
		break;
	};
}
