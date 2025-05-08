# Define version or set it from git

IF (NOT YADE_VERSION)
  IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/RELEASE )
    #Release file is found
    SET(READFILE cat)
    execute_process(
      COMMAND ${READFILE} "RELEASE"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE YADE_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  ELSEIF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    #Use git for version defining
    execute_process(
      COMMAND git log -n1 --pretty=oneline
      COMMAND cut -c1-7
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE VERSION_GIT
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
      COMMAND git log -n1 --pretty=fuller --date=iso
      COMMAND grep AuthorDate
      COMMAND cut -c13-22
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE VERSION_DATE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET(YADE_VERSION "${VERSION_DATE}.git-${VERSION_GIT}")

    #git log -n1 --pretty=format:"%ai_%h"
  ELSEIF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.bzr/branch/last-revision)
    #Use bzr for version defining
    execute_process(
      COMMAND less last-revision
      COMMAND cut -c13-20
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.bzr/branch/
      OUTPUT_VARIABLE VERSION_GIT
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
      COMMAND bzr log -l 1 --gnu-changelog
      COMMAND head -n 1
      COMMAND cut -c1-10
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE VERSION_DATE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET(YADE_VERSION "${VERSION_DATE}.git-${VERSION_GIT}")
  ELSE (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git )
    SET (YADE_VERSION "Unknown")
  ENDIF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/RELEASE )
ENDIF (NOT YADE_VERSION)


MESSAGE (STATUS "Version is set to " ${YADE_VERSION})
