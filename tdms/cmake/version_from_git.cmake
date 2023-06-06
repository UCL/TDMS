find_package(Git)
if(GIT_FOUND)

	# Check we're in a git repo.
	execute_process(OUTPUT_VARIABLE IS_REPO
		COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		ERROR_QUIET)

	if (NOT IS_REPO)
		message(WARNING
		       	"No git repository found. You can safely ignore this warning if you're compiling from a tar of the source code.")
		unset(TDMS_VERSION)
		return()
	endif()

	# Get current branch name.
	execute_process(OUTPUT_VARIABLE CURRENT_BRANCH
		COMMAND ${GIT_EXECUTABLE} symbolic-ref --short HEAD
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_STRIP_TRAILING_WHITESPACE )

	# Get current tag (or short-form of the commit hash if not tagged).
	execute_process(OUTPUT_VARIABLE TAG_OR_SHORT_HASH
		COMMAND ${GIT_EXECUTABLE} describe --tags --always
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_STRIP_TRAILING_WHITESPACE )

	string(REGEX MATCH "v[0-9]+\.[0-9]+\.[0-9]+" IS_SEMVER "${TAG_OR_SHORT_HASH}")
	if (CURRENT_BRANCH STREQUAL "main" OR CURRENT_BRANCH STREQUAL "" AND IS_SEMVER)
	    # If we're on main or in a detatched head, and we're tagged with a
	    # sememantic version (of the form v1.2.3) then that is the version.
		set(TDMS_VERSION "${TAG_OR_SHORT_HASH}")
	else()
		# If it's not a tag on main or if it's a commit hash then the "version"
		# placeholder is <branchname>_<tag>.
		#
		# E.g. my-feature-branch_04ad23
		#      my-feature-branch_taggingForSomeReason
		set(TDMS_VERSION "${CURRENT_BRANCH}_${TAG_OR_SHORT_HASH}")
	endif()

endif()
