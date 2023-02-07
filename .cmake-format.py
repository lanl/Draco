# -------------------------------------------*-python-*------------------------------------------- #
# \file   .cmake-format.py
# \author Kelly Thompson <kgt@lanl.gov>
# \brief  Options affecting listfile parsing
# \note   Copyright (C) 2021-2023 Triad National Security, LLC., All rights reserved.
#
# https://cmake-format.readthedocs.io/en/latest/configuration.html
# ------------------------------------------------------------------------------------------------ #
with section("parse"):  # noqa: F821

    # Specify structure for custom cmake functions
    additional_commands = {

        "add_component_library": {
            "pargs": 1,
            "flags": ["NOEXPORT", "PROVIDE_DLL_DEPS"],
            "kwargs": {"TARGET": '*', "LIBRARY_NAME": '*', "SOURCES": '*', "HEADERS": '*',
                       "INCLUDE_DIRS": '*', "TARGET_DEPS": '*', "LANGUAGE": '*',
                       "LINK_LANGUAGE": '*', "LIBRARY_TYPE": '*'}},

        "add_component_executable": {
            "kwargs": {"TARGET": '*', "LIBRARY_NAME": '*', "SOURCES": '*', "HEADERS": '*',
                       "INCLUDE_DIRS": '*', "TARGET_DEPS": '*', "EXE_NAME": '*', "FOLDER": '*'}},

        "add_scalar_tests": {
            "kwargs": {"SOURCES": '*', "DEPS": '*', "TEST_ARGS": '*', "PASS_REGEX": '*', "LABELS":
                       '*'}},

        "add_parallel_tests": {
            "flags": ["MPI_PLUS_OMP"],
            "kwargs": {"SOURCES": '*', "DEPS": '*', "TEST_ARGS": '*', "PASS_REGEX": '*',
                       "FAIL_REGEX": '*', "LABEL": '*', "PE_LIST": '*'}},

        "add_app_unit_test": {
            "kwargs": {"DRIVER": '*', "APPS": '*', "TEST_ARGS": '*', "LABELS": '*', "APP": '*',
                       "PE_LIST": '*', "STDINFILE": '*', "GOLDFILE": '*'}},

        "cmake_add_fortran_subdirectory": {
            "pargs": 1,
            "flags": ["VERBOSE", "NO_EXTERNAL_INSTALL"],
            "kwargs": {"PROJECT": '*', "ARCHIVE_DIR": '*', "RUNTIME_DIR": '*', "LIBRARIES": '*',
                       "TARGET_NAMES": '*', "DEPENDS": '*', "CMAKE_COMMAND_LINE": '*'}},

        "cmake_parse_arguments": {
            "pargs": 5,
            "flags": ["VERBOSE"],
            "kwargs": {"FOOBAR": '*'}},

        "fortrancinterface_header": {
            "pargs": 1,
            "kwargs": {"MACRO_NAMESPACE": '*', "SYMBOLS": '*', "SYMBOL_NAMESPACE": '*'}},

        "register_parallel_test": {
            "pargs": 0,
            "flags": ["FOO_PROPERTIES", "NO_EXTERNAL_INSTALL"],
            "kwargs": {"CMD_ARGS": '*', "COMMAND": '*', "NUMPE": '*', "TARGET": '*'}},

        "register_scalar_test": {
            "pargs": 0,
            "flags": ["FOO_PROPERTIES", "NO_EXTERNAL_INSTALL"],
            "kwargs": {"TARGET": '*', "COMMAND": '*', "CMD_ARGS": '*'}},

        "set_target_properties": {
            "pargs": 1,
            "flags": ["FOO_PROPERTIES", "NO_EXTERNAL_INSTALL"],
            "kwargs": {"COMPILE_DEFINITIONS": '*', "DEPENDS": '*', "FOLDER": '*',
                       "IMPORTED_CONFIGURATIONS": '*',
                       "IMPORTED_IMPLIB": '*',
                       "IMPORTED_IMPLIB_DEBUG": '*',
                       "IMPORTED_IMPLIB_RELEASE": '*',
                       "IMPORTED_LINK_INTERFACE_LANGUAGES": '*',
                       "IMPORTED_LOCATION_DEBUG": '*',
                       "IMPORTED_LOCATION_RELEASE": '*', "IMPORTED_LINK_INTERFACE_LIBRARIES": '*',
                       "IMPORTED_LINK_INTERFACE_LANGUAGES": '*',
                       "IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG": '*',
                       "IMPORTED_LOCATION": '*',
                       "INTERFACE_INCLUDE_DIRECTORIES": '*',
                       "LINKER_LANGUAGE": '*',
                       "OUTPUT_NAME": '*',
                       "RUNTIME_OUTPUT_DIRECTORY": '*',
                       "VS_KEYWORD": '*',
                       "WINDOWS_EXPORT_ALL_SYMBOLS": '*',
                       "XCODE_ATTRIBUTE_ENABLE_OPENMP_SUPPORT": '*'}},

        "set_tests_properties": {
            "pargs": 1,
            "flags": ["FOO_PROPERTIES", "NO_EXTERNAL_INSTALL"],
            "kwargs": {"DEPENDS": '*', "ENVIRONMENT": '*', "FAIL_REGULAR_EXPRESSION": '*',
                       "LABELS": '*', "PASS_REGULAR_EXPRESSION": '*', "PROCESSORS": '*',
                       "RESOURCE_LOCK": '*', "WORKING_DIRECTORY": '*'}},
    }

    # Specify variable tags.
    vartags = []

    # Specify property tags.
    proptags = []

# -----------------------------
# Options affecting formatting.
# -----------------------------
with section("format"):  # noqa: F821

    # How wide to allow formatted cmake files
    line_width = 100

    # How many spaces to tab for indent
    tab_size = 2

    # If an argument group contains more than this many sub-groups (parg or kwarg groups) then force
    # it to a vertical layout.
    max_subgroups_hwrap = 2

    # If a positional argument group contains more than this many arguments, then force it to a
    # vertical layout.
    max_pargs_hwrap = 6

    # If a cmdline positional group consumes more than this many lines without nesting, then
    # invalidate the layout (and nest)
    max_rows_cmdline = 2

    # If true, separate flow control names from their parentheses with a space
    separate_ctrl_name_with_space = False

    # If true, separate function names from parentheses with a space
    separate_fn_name_with_space = False

    # If a statement is wrapped to more than one line, than dangle the closing parenthesis on its
    # own line.
    dangle_parens = False

    # If the trailing parenthesis must be 'dangled' on its on line, then align it to this reference:
    # `prefix`: the start of the statement, `prefix-indent`: the start of the statement, plus one
    # indentation level, `child`: align to the column of the arguments
    dangle_align = 'prefix'

    # If the statement spelling length (including space and parenthesis) is smaller than this
    # amount, then force reject nested layouts.
    min_prefix_chars = 4

    # If the statement spelling length (including space and parenthesis) is larger than the tab
    # width by more than this amount, then force reject un-nested layouts.
    max_prefix_chars = 10

    # If a candidate layout is wrapped horizontally but it exceeds this many lines, then reject the
    # layout.
    max_lines_hwrap = 2

    # What style line endings to use in the output.
    line_ending = 'unix'

    # Format command names consistently as 'lower' or 'upper' case
    command_case = 'canonical'

    # Format keywords consistently as 'lower' or 'upper' case
    keyword_case = 'unchanged'

    # A list of command names which should always be wrapped
    always_wrap = []

    # If true, the argument lists which are known to be sortable will be sorted lexicographicall
    enable_sort = True

    # If true, the parsers may infer whether or not an argument list is sortable (without
    # annotation).
    autosort = True

    # By default, if cmake-format cannot successfully fit everything into the desired linewidth it
    # will apply the last, most agressive attempt that it made. If this flag is True, however,
    # cmake-format will print error, exit with non-zero status code, and write-out nothing
    require_valid_layout = False

    # A dictionary mapping layout nodes to a list of wrap decisions. See the documentation for more
    # information.
    layout_passes = {}

# ------------------------------------------------
# Options affecting comment reflow and formatting.
# ------------------------------------------------
with section("markup"):  # noqa: F821

    # What character to use for bulleted lists
    bullet_char = '*'

    # What character to use as punctuation after numerals in an enumerated list
    enum_char = '.'

    # If comment markup is enabled, don't reflow the first comment block in each listfile. Use this
    # to preserve formatting of your copyright/license statements.
    first_comment_is_literal = True

    # If comment markup is enabled, don't reflow any comment block which matches this (regex)
    # pattern. Default is `None` (disabled).
    literal_comment_pattern = None

    # Regular expression to match preformat fences in comments default=
    # ``r'^\s*([`~]{3}[`~]*)(.*)$'``
    fence_pattern = '^\\s*([`~]{3}[`~]*)(.*)$'

    # Regular expression to match rulers in comments default= ``r'^\s*[^\w\s]{3}.*[^\w\s]{3}$'``
    ruler_pattern = '^\\s*[^\\w\\s]{3}.*[^\\w\\s]{3}$'

    # If a comment line matches starts with this pattern then it is explicitly a trailing comment
    # for the preceeding argument. Default is '#<'
    explicit_trailing_pattern = '#<'

    # If a comment line starts with at least this many consecutive hash characters, then don't
    # lstrip() them off. This allows for lazy hash rulers where the first hash char is not separated
    # by space
    hashruler_min_length = 10

    # If true, then insert a space between the first hash char and remaining hash chars in a hash
    # ruler, and normalize its length to fill the column
    canonicalize_hashrulers = True

    # enable comment markup parsing and reflow
    enable_markup = True

# ----------------------------
# Options affecting the linter
# ----------------------------
with section("lint"):  # noqa: F821

    # a list of lint codes to disable
    disabled_codes = ["C0103"]

    # regular expression pattern describing valid function names
    function_pattern = '[0-9a-z_]+'

    # regular expression pattern describing valid macro names
    macro_pattern = '[0-9A-Z_]+'

    # regular expression pattern describing valid names for variables with global (cache) scope
    global_var_pattern = '[A-Z][0-9A-Z_]+'

    # regular expression pattern describing valid names for variables with global scope (but
    # internal semantic)
    internal_var_pattern = '_[A-Z][0-9A-Z_]+'

    # regular expression pattern describing valid names for variables with local scope
    local_var_pattern = '[a-z][a-z0-9_]+'

    # regular expression pattern describing valid names for privatedirectory variables
    private_var_pattern = '_[0-9a-z_]+'

    # regular expression pattern describing valid names for public directory variables
    public_var_pattern = '[A-Z][0-9A-Z_]+'

    # regular expression pattern describing valid names for function/macro arguments and loop
    # variables.
    argument_var_pattern = '[a-z][a-z0-9_]+'

    # regular expression pattern describing valid names for keywords used in functions or macros
    keyword_pattern = '[A-Z][0-9A-Z_]+'

    # In the heuristic for C0201, how many conditionals to match within a loop in before considering
    # the loop a parser.
    max_conditionals_custom_parser = 2

    # Require at least this many newlines between statements
    min_statement_spacing = 1

    # Require no more than this many newlines between statements
    max_statement_spacing = 2
    max_returns = 6
    max_branches = 12
    max_arguments = 5
    max_localvars = 15
    max_statements = 50

# -------------------------------
# Options affecting file encoding
# -------------------------------
with section("encode"):  # noqa: F821

    # If true, emit the unicode byte-order mark (BOM) at the start of the file
    emit_byteorder_mark = False

    # Specify the encoding of the input file. Defaults to utf-8
    input_encoding = 'utf-8'

    # Specify the encoding of the output file. Defaults to utf-8. Note that cmake only claims to
    # support utf-8 so be careful when using anything else
    output_encoding = 'utf-8'

# -------------------------------------
# Miscellaneous configurations options.
# -------------------------------------
with section("misc"):  # noqa: F821

    # A dictionary containing any per-command configuration overrides. Currently only `command_case`
    # is supported.
    per_command = {}
