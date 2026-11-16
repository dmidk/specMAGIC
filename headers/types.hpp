#pragma once



// Type used for BOTH output radiation products 
// AND for storing the incoming satellite image 
// The user MUST verify that this is a valid type for satellite image
// before changing.
using MAGIC_INT = short int;

// The data type used for calculation
#if defined(MAGIC_PRECISION_HI)
using MAGIC_REAL = double;
#elif defined(MAGIC_PRECISION_LO)
using MAGIC_REAL = float;
#else
#error "Precision not specified!"
#endif

// For things that must ALWAYS be in high precision
using MAGIC_EXACT = double;

// Whether or not a given pixel should be processed
enum class PixelState { Night = 0, Day = 1, Invalid = 99 };

// To distinguish between GHI and DNI
enum class IrradianceMode { Beam, Global };

// The source of the albedo data being used
// INPUT is an albedo field given for the time of the image (InputFields), such as LSA SAF's
enum class AlbedoType { LANDMAP, MODIS, INPUT };
