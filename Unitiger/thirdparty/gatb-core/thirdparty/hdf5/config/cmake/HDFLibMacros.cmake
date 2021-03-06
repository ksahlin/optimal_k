#-------------------------------------------------------------------------------
MACRO (EXTERNAL_JPEG_LIBRARY compress_type libtype jpeg_pic)
  # May need to build JPEG with PIC on x64 machines with gcc
  # Need to use CMAKE_ANSI_CFLAGS define so that compiler test works

  IF (${compress_type} MATCHES "SVN")
    EXTERNALPROJECT_ADD (JPEG
        SVN_REPOSITORY ${JPEG_URL}
        # [SVN_REVISION rev] 
        INSTALL_COMMAND ""
        CMAKE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF_PACKAGE_EXT:STRING=${HDF_PACKAGE_EXT}
            -DJPEG_EXTERNALLY_CONFIGURED:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            -DCMAKE_ANSI_CFLAGS:STRING=${jpeg_pic}
    ) 
  ELSEIF (${compress_type} MATCHES "TGZ")
    EXTERNALPROJECT_ADD (JPEG
        URL ${JPEG_URL}
        URL_MD5 ""
        INSTALL_COMMAND ""
        CMAKE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF_PACKAGE_EXT:STRING=${HDF_PACKAGE_EXT}
            -DJPEG_EXTERNALLY_CONFIGURED:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            -DCMAKE_ANSI_CFLAGS:STRING=${jpeg_pic}
    ) 
  ENDIF (${compress_type} MATCHES "SVN")
  EXTERNALPROJECT_GET_PROPERTY (JPEG BINARY_DIR SOURCE_DIR) 

  IF (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    IF (WIN32)
      SET (JPEG_LIB_NAME "jpeg_D")
    ELSE (WIN32)
      SET (JPEG_LIB_NAME "jpeg_debug")
    ENDIF (WIN32)
  ELSE (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET (JPEG_LIB_NAME "jpeg")
  ENDIF (${CMAKE_BUILD_TYPE} MATCHES "Debug")

  # Create imported target szip
  ADD_LIBRARY(jpeg ${libtype} IMPORTED)
  ADD_DEPENDENCIES (jpeg JPEG)

  IF (${libtype} MATCHES "SHARED")
    IF (WIN32)
      IF (MINGW)
        SET_TARGET_PROPERTIES(jpeg PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${JPEG_LIB_NAME}.lib"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${JPEG_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ELSE (MINGW)
        SET_TARGET_PROPERTIES(jpeg PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${CMAKE_IMPORT_LIBRARY_PREFIX}${JPEG_LIB_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${CMAKE_IMPORT_LIBRARY_PREFIX}${JPEG_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ENDIF (MINGW)
    ELSE (WIN32)
      IF (CYGWIN)
        SET_TARGET_PROPERTIES(jpeg PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_IMPORT_LIBRARY_PREFIX}${JPEG_LIB_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_IMPORT_LIBRARY_PREFIX}${JPEG_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ELSE (CYGWIN)
        SET_TARGET_PROPERTIES(jpeg PROPERTIES
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}${JPEG_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            IMPORTED_SONAME "${CMAKE_SHARED_LIBRARY_PREFIX}${JPEG_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}.${JPEG_VERSION_STRING}"
            SOVERSION "${JPEG_VERSION_STRING}"
        )
      ENDIF (CYGWIN)
    ENDIF (WIN32)
  ELSE (${libtype} MATCHES "SHARED")
    IF (WIN32 AND NOT MINGW)
      SET_TARGET_PROPERTIES(jpeg PROPERTIES
          IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/lib${JPEG_LIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      )
    ELSE (WIN32 AND NOT MINGW)
      SET_TARGET_PROPERTIES(jpeg PROPERTIES
          IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${JPEG_LIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      )
    ENDIF (WIN32 AND NOT MINGW)
  ENDIF (${libtype} MATCHES "SHARED")

#  INCLUDE (${BINARY_DIR}/JPEG-targets.cmake)  
  SET (JPEG_LIBRARY "jpeg")
  
  SET (JPEG_INCLUDE_DIR_GEN "${BINARY_DIR}")
  SET (JPEG_INCLUDE_DIR "${SOURCE_DIR}/src")
  SET (JPEG_FOUND 1)
  SET (JPEG_LIBRARIES ${JPEG_LIBRARY})
  SET (JPEG_INCLUDE_DIRS ${JPEG_INCLUDE_DIR_GEN} ${JPEG_INCLUDE_DIR})
ENDMACRO (EXTERNAL_JPEG_LIBRARY)

#-------------------------------------------------------------------------------
MACRO (PACKAGE_JPEG_LIBRARY compress_type)
  ADD_CUSTOM_TARGET (JPEG-GenHeader-Copy ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${JPEG_INCLUDE_DIR_GEN}/jconfig.h ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
      COMMENT "Copying ${JPEG_INCLUDE_DIR_GEN}/jconfig.h to ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/"
  )
  SET (EXTERNAL_HEADER_LIST ${EXTERNAL_HEADER_LIST} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/jconfig.h)
  IF (${compress_type} MATCHES "SVN" OR ${compress_type} MATCHES "TGZ")
    ADD_DEPENDENCIES (JPEG-GenHeader-Copy JPEG)
  ENDIF (${compress_type} MATCHES "SVN" OR ${compress_type} MATCHES "TGZ")
ENDMACRO (PACKAGE_JPEG_LIBRARY)

#-------------------------------------------------------------------------------
MACRO (EXTERNAL_SZIP_LIBRARY compress_type libtype encoding)
  IF (${compress_type} MATCHES "SVN")
    EXTERNALPROJECT_ADD (SZIP
        SVN_REPOSITORY ${SZIP_URL}
        # [SVN_REVISION rev] 
        INSTALL_COMMAND ""
        CMAKE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF_PACKAGE_EXT:STRING=${HDF_PACKAGE_EXT}
            -DSZIP_EXTERNALLY_CONFIGURED:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            -DCMAKE_ANSI_CFLAGS:STRING=${CMAKE_ANSI_CFLAGS}
            -DSZIP_ENABLE_ENCODING:BOOL=${encoding}
    ) 
  ELSEIF (${compress_type} MATCHES "TGZ")
    EXTERNALPROJECT_ADD (SZIP
        URL ${SZIP_URL}
        URL_MD5 ""
        INSTALL_COMMAND ""
        CMAKE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF_PACKAGE_EXT:STRING=${HDF_PACKAGE_EXT}
            -DSZIP_EXTERNALLY_CONFIGURED:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            -DCMAKE_ANSI_CFLAGS:STRING=${CMAKE_ANSI_CFLAGS}
            -DSZIP_ENABLE_ENCODING:BOOL=${encoding}
    ) 
  ENDIF (${compress_type} MATCHES "SVN")
  EXTERNALPROJECT_GET_PROPERTY (SZIP BINARY_DIR SOURCE_DIR) 

  IF (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    IF (WIN32)
      SET (SZIP_LIB_NAME "szip_D")
    ELSE (WIN32)
      SET (SZIP_LIB_NAME "szip_debug")
    ENDIF (WIN32)
  ELSE (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET (SZIP_LIB_NAME "szip")
  ENDIF (${CMAKE_BUILD_TYPE} MATCHES "Debug")

  # Create imported target szip
  ADD_LIBRARY(szip ${libtype} IMPORTED)
  ADD_DEPENDENCIES (szip SZIP)

  IF (${libtype} MATCHES "SHARED")
    IF (WIN32)
      IF (MINGW)
        SET_TARGET_PROPERTIES(szip PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${SZIP_LIB_NAME}.lib"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${SZIP_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ELSE (MINGW)
        SET_TARGET_PROPERTIES(szip PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${CMAKE_IMPORT_LIBRARY_PREFIX}${SZIP_LIB_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${CMAKE_IMPORT_LIBRARY_PREFIX}${SZIP_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ENDIF (MINGW)
    ELSE (WIN32)
      IF (CYGWIN)
        SET_TARGET_PROPERTIES(szip PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_IMPORT_LIBRARY_PREFIX}${SZIP_LIB_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_IMPORT_LIBRARY_PREFIX}${SZIP_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ELSE (CYGWIN)
        SET_TARGET_PROPERTIES(szip PROPERTIES
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}${SZIP_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            IMPORTED_SONAME "${CMAKE_SHARED_LIBRARY_PREFIX}${SZIP_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}.${SZIP_VERSION_STRING}"
            SOVERSION "${SZIP_VERSION_STRING}"
        )
      ENDIF (CYGWIN)
    ENDIF (WIN32)
  ELSE (${libtype} MATCHES "SHARED")
    IF (WIN32 AND NOT MINGW)
      SET_TARGET_PROPERTIES(szip PROPERTIES
          IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/lib${SZIP_LIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      )
    ELSE (WIN32 AND NOT MINGW)
      SET_TARGET_PROPERTIES(szip PROPERTIES
          IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${SZIP_LIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      )
    ENDIF (WIN32 AND NOT MINGW)
  ENDIF (${libtype} MATCHES "SHARED")

#  INCLUDE (${BINARY_DIR}/SZIP-targets.cmake)  
  SET (SZIP_LIBRARY "szip")

  SET (SZIP_INCLUDE_DIR_GEN "${BINARY_DIR}")
  SET (SZIP_INCLUDE_DIR "${SOURCE_DIR}/src")
  SET (SZIP_FOUND 1)
  SET (SZIP_LIBRARIES ${SZIP_LIBRARY})
  SET (SZIP_INCLUDE_DIRS ${SZIP_INCLUDE_DIR_GEN} ${SZIP_INCLUDE_DIR})
ENDMACRO (EXTERNAL_SZIP_LIBRARY)

#-------------------------------------------------------------------------------
MACRO (PACKAGE_SZIP_LIBRARY compress_type)
  ADD_CUSTOM_TARGET (SZIP-GenHeader-Copy ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SZIP_INCLUDE_DIR_GEN}/SZconfig.h ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
      COMMENT "Copying ${SZIP_INCLUDE_DIR_GEN}/SZconfig.h to ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/"
  )
  SET (EXTERNAL_HEADER_LIST ${EXTERNAL_HEADER_LIST} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/SZconfig.h)
  IF (${compress_type} MATCHES "SVN" OR ${compress_type} MATCHES "TGZ")
    ADD_DEPENDENCIES (SZIP-GenHeader-Copy SZIP)
  ENDIF (${compress_type} MATCHES "SVN" OR ${compress_type} MATCHES "TGZ")
ENDMACRO (PACKAGE_SZIP_LIBRARY)

#-------------------------------------------------------------------------------
MACRO (EXTERNAL_ZLIB_LIBRARY compress_type libtype)
  IF (${compress_type} MATCHES "SVN")
    EXTERNALPROJECT_ADD (ZLIB
        SVN_REPOSITORY ${ZLIB_URL}
        # [SVN_REVISION rev] 
        INSTALL_COMMAND ""
        CMAKE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF_PACKAGE_EXT:STRING=${HDF_PACKAGE_EXT}
            -DZLIB_EXTERNALLY_CONFIGURED:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            -DCMAKE_ANSI_CFLAGS:STRING=${CMAKE_ANSI_CFLAGS}
    ) 
  ELSEIF (${compress_type} MATCHES "TGZ")
    EXTERNALPROJECT_ADD (ZLIB
        URL ${ZLIB_URL}
        URL_MD5 ""
        INSTALL_COMMAND ""
        CMAKE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF_PACKAGE_EXT:STRING=${HDF_PACKAGE_EXT}
            -DZLIB_EXTERNALLY_CONFIGURED:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            -DCMAKE_ANSI_CFLAGS:STRING=${CMAKE_ANSI_CFLAGS}
    ) 
  ENDIF (${compress_type} MATCHES "SVN")
  EXTERNALPROJECT_GET_PROPERTY (ZLIB BINARY_DIR SOURCE_DIR) 

  IF (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    IF (WIN32)
      SET (ZLIB_LIB_NAME "zlib_D")
    ELSE (WIN32)
      SET (ZLIB_LIB_NAME "z_debug")
    ENDIF (WIN32)
  ELSE (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    IF (WIN32)
      SET (ZLIB_LIB_NAME "zlib")
    ELSE (WIN32)
      SET (ZLIB_LIB_NAME "z")
    ENDIF (WIN32)
  ENDIF (${CMAKE_BUILD_TYPE} MATCHES "Debug")

  # Create imported target szip
  ADD_LIBRARY(zlib ${libtype} IMPORTED)
  ADD_DEPENDENCIES (zlib ZLIB)
  
  IF (${libtype} MATCHES "SHARED")
    IF (WIN32)
      IF (MINGW)
        SET_TARGET_PROPERTIES(zlib PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${ZLIB_LIB_NAME}.lib"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${ZLIB_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ELSE (MINGW)
        SET_TARGET_PROPERTIES(zlib PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${CMAKE_IMPORT_LIBRARY_PREFIX}${ZLIB_LIB_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${CMAKE_IMPORT_LIBRARY_PREFIX}${ZLIB_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ENDIF (MINGW)
    ELSE (WIN32)
      IF (CYGWIN)
        SET_TARGET_PROPERTIES(zlib PROPERTIES
            IMPORTED_IMPLIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_IMPORT_LIBRARY_PREFIX}${ZLIB_LIB_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}"
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_IMPORT_LIBRARY_PREFIX}${ZLIB_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        )
      ELSE (CYGWIN)
        SET_TARGET_PROPERTIES(zlib PROPERTIES
            IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}${ZLIB_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            IMPORTED_SONAME "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}${ZLIB_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}.${ZLIB_VERSION_STRING}"
            SOVERSION "${ZLIB_VERSION_STRING}"
        )
      ENDIF (CYGWIN)
    ENDIF (WIN32)
  ELSE (${libtype} MATCHES "SHARED")
    IF (WIN32 AND NOT MINGW)
      SET_TARGET_PROPERTIES(zlib PROPERTIES
          IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/lib${ZLIB_LIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      )
    ELSE (WIN32 AND NOT MINGW)
      SET_TARGET_PROPERTIES(zlib PROPERTIES
          IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${ZLIB_LIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      )
    ENDIF (WIN32 AND NOT MINGW)
  ENDIF (${libtype} MATCHES "SHARED")

#  INCLUDE (${BINARY_DIR}/ZLIB-targets.cmake)  
  SET (ZLIB_LIBRARY "zlib")
  
  SET (ZLIB_INCLUDE_DIR_GEN "${BINARY_DIR}")
  SET (ZLIB_INCLUDE_DIR "${SOURCE_DIR}")
  SET (ZLIB_FOUND 1)
  SET (ZLIB_LIBRARIES ${ZLIB_LIBRARY})
  SET (ZLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR_GEN} ${ZLIB_INCLUDE_DIR})
ENDMACRO (EXTERNAL_ZLIB_LIBRARY)

#-------------------------------------------------------------------------------
MACRO (PACKAGE_ZLIB_LIBRARY compress_type)
  ADD_CUSTOM_TARGET (ZLIB-GenHeader-Copy ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${ZLIB_INCLUDE_DIR_GEN}/zconf.h ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
      COMMENT "Copying ${ZLIB_INCLUDE_DIR_GEN}/zconf.h to ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/"
  )
  SET (EXTERNAL_HEADER_LIST ${EXTERNAL_HEADER_LIST} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/zconf.h)
  IF (${compress_type} MATCHES "SVN" OR ${compress_type} MATCHES "TGZ")
    ADD_DEPENDENCIES (ZLIB-GenHeader-Copy ZLIB)
  ENDIF (${compress_type} MATCHES "SVN" OR ${compress_type} MATCHES "TGZ")
ENDMACRO (PACKAGE_ZLIB_LIBRARY)
