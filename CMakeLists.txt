cmake_minimum_required(VERSION 3.5)

project(SlicerBoneRegistration)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/VertRegExtension")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Lance Levine")
set(EXTENSION_DESCRIPTION "This extension will generate a DRR and register that to a fluoro image.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/VertRegExtension.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/VertRegExtension/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(VertReg)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
