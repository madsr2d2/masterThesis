"""
Shared data-loading functions for the benzyl alcohol / H2O2 kinetics dataset.

The extraction functions below (parse_experiment_data through
find_substrate_type, plus populate_experimental_data_from_directory) are
copied verbatim from masterThesis.ipynb (cells 3, 5-12) so there is exactly
one implementation of the parsing logic instead of one per script. If the
notebook's own copies ever change, update both -- see DATA_VERIFICATION.md
for why keeping this logic byte-identical to what actually produced
data/experiment_data.csv matters.

load_experiment() is new: a per-experiment convenience wrapper (metadata +
raw time series together) for tools like plot_kinetics.py that don't need
the full-dataset DataFrame.
"""
import re
from pathlib import Path
import pandas as pd
import numpy as np
import os

SUBSTRATE_PROPERTIES = {
    "BnOH": {"abs": 285, "e": 1.23},
    "4OMe-BnOH": {"abs": 300, "e": 7.53},
}


# ---------------------------------------------------------------------------
# Per-experiment extraction corrections
# ---------------------------------------------------------------------------
# Experiments 32, 34, 35, 36 and 37 are buffer-concentration titrations whose
# .xls sheets lay out EIGHT planned cuvettes -- rows 1-4 with enzyme, rows 5-8
# without -- while only four channels were ever measured. These were no-enzyme
# days, so the four cuvettes that actually ran were 5-8; but
# find_numeric_values_below_header reads the first four rows of the table and
# therefore picks up the with-enzyme plan rows. Two columns come out wrong:
#
#   [enz]  extracted as 0.241 / 0.270 mM. The runs were enzyme-free -- the
#          filenames say "with_NO_E" and all five sit in the hand-sorted
#          data/Mads/"No enzyme"/ folder. 22 of the 27 experiments in that
#          folder extract correctly as 0; these five are the exceptions.
#   [buf]  extracted as a flat 50 mM, because every cuvette receives the same
#          buffer VOLUME (1 ml into 2 ml total) and only the stock differs --
#          and the stock appears solely as a text label in the "kuv" column
#          ("1 (0.1M)", "2 (0.2M)", ...), which the volume-based extraction
#          cannot see. The cuvette concentration is half the stock.
#
# Corrected, these five become enzyme-free buffer titrations spanning
# 3.125-200 mM at constant substrate. See DATA_VERIFICATION.md, 2026-08-30.
EXPERIMENT_CORRECTIONS = {
    # experiment: {column: scalar applied to every sample, or per-sample list}
    32: {"[enz]": 0.0, "[buf]": [50.0, 100.0, 150.0, 200.0]},
    34: {"[enz]": 0.0, "[buf]": [25.0, 12.5, 6.25, 3.125]},
    35: {"[enz]": 0.0, "[buf]": [50.0, 100.0, 150.0, 200.0]},
    36: {"[enz]": 0.0, "[buf]": [50.0, 100.0, 150.0, 200.0]},
    37: {"[enz]": 0.0, "[buf]": [50.0, 100.0, 150.0, 200.0]},

    # Exps 79 and 80 are enzyme runs whose [enz] extracted as zero. Their
    # cuvette tables DO carry an "[Enz] mmol/l" column, but every measured row
    # holds 0.000001 -- a broken formula, about 14,000x too low -- which the
    # extraction read faithfully and rounding to 3 dp turned into 0.0.
    #
    # The right value is in the sheet's header block, on the "kuv" row of the
    # enzyme stock calculation, and it checks out against the volumes: the stock
    # is 0.559618 mM, and 0.559618 * 0.05/2 = 0.01399 for exp 79, * 0.1/2 =
    # 0.027981 for exp 80. Both filenames say "with_E". (That header row is the
    # same one that matches the table exactly in every healthy sheet -- exp 2
    # declares kuv = 0.17533 and its table column reads 0.17533.)
    #
    # Checked across all 98: 63 sheets declare a header kuv, 58 agree with the
    # compiled [enz], and the only other disagreements are exps 32 and 34-37,
    # where the kuv belongs to the planned-but-unmeasured with-enzyme rows. So
    # the broken column is confined to these two.
    #
    # Exp 80 is IN USE, so until now an enzyme run was sitting in the dataset
    # indistinguishable from an enzyme-free control -- precisely the set the
    # catalyst-independent rate constants are meant to be fitted on.
    # Ruled 2026-08-30; see DATA_VERIFICATION.md.
    79: {"[enz]": 0.014},
    80: {"[enz]": 0.028},
}


def apply_experiment_corrections(experiment_number, sample_index, concentrations):
    """
    Overrides extracted concentrations for the experiments whose .xls layout
    defeats the generic extraction (see EXPERIMENT_CORRECTIONS above).

    Args:
        experiment_number (int): The experiment number.
        sample_index (int): Zero-based index of the sample within the experiment.
        concentrations (dict): Extracted {"[enz]", "[buf]", "[h2o2]", "[sub]"}.

    Returns:
        dict: A new dict with any corrections applied.
    """
    corrected = dict(concentrations)
    for column, value in (EXPERIMENT_CORRECTIONS.get(experiment_number) or {}).items():
        if isinstance(value, (list, tuple)):
            if sample_index < len(value):
                corrected[column] = value[sample_index]
        else:
            corrected[column] = value
    return corrected


def parse_experiment_data(file_path):
    """
    Parses the experiment data file and extracts the experiment number, date,
    and time series with values for each sample.

    Args:
        file_path (str): Path to the experiment data file.

    Returns:
        dict: Parsed data in the format:
            {
                "num": experiment_number,
                "date": date_collected,
                "samples": {
                    "Sample001": {"time": [times], "values": [values]},
                    ...
                }
            }
    """
    parsed_data = {}

    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Extract the experiment number
        experiment_line = lines[0].strip()
        match = re.search(r"(?:rate|mads_t)(\d+)\.rre", experiment_line, re.IGNORECASE)
        if match:
            parsed_data["num"] = int(match.group(1))
        else:
            parsed_data["num"] = None
            print("Warning: Experiment number not found.")

        # Extract the date
        if len(lines) > 5:
            date_line = lines[5].strip().split('\t')
            if len(date_line) > 1:
                parsed_data["date"] = date_line[1]
            else:
                parsed_data["date"] = None
                print("Warning: Date not found in the expected format.")
        else:
            parsed_data["date"] = None
            print("Warning: Date line is missing.")

        # Extract sample names (case-insensitive normalization)
        if len(lines) > 2:
            sample_names = lines[2].strip().split('\t')[1::2]
            normalized_sample_names = {sample.lower(): sample for sample in sample_names}
            parsed_data["samples"] = {normalized_sample_names[sample]: {"time": [], "values": []}
                                      for sample in normalized_sample_names}
        else:
            parsed_data["samples"] = {}
            print("Warning: Sample names line is missing.")

        # Locate the start of the data section
        for i, line in enumerate(lines):
            if line.lower().startswith("ss.sss"):
                data_start_line = i + 1
                break
        else:
            raise ValueError("Data section not found in the file.")

        # Parse the data
        for line in lines[data_start_line:]:
            values = line.strip().split('\t')
            if not values or len(values) < 2:  # Skip empty or incomplete lines
                continue
            try:
                time = float(values[0])  # Extract time from the first column
            except ValueError:
                print(f"Warning: Invalid time value in line: {line}")
                continue

            for i, sample in enumerate(normalized_sample_names):
                value_index = 2 * i + 1  # Locate the corresponding value column
                if value_index < len(values):  # Check if the value exists
                    try:
                        value = float(values[value_index])
                        parsed_data["samples"][normalized_sample_names[sample]]["time"].append(time)
                        parsed_data["samples"][normalized_sample_names[sample]]["values"].append(value)
                    except ValueError:
                        print(f"Warning: Invalid value for {sample} in line: {line}")
                        continue

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None
    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")
        return None

    return parsed_data


def find_and_parse_experiment_file(experiment_number, directory, sheet_name=None):
    """
    Finds the associated .xls file for a given experiment number by searching recursively
    through subdirectories starting from the provided path and parses the specified sheet.
    If the primary method fails, uses the #num pattern as a fallback.

    Args:
        experiment_number (int): The experiment number to search for (e.g., 26 for mads_t026).
        directory (str): The directory where the search starts.
        sheet_name (str, optional): Name of the sheet to parse. If None, the first sheet is used.

    Returns:
        tuple: A tuple containing the filename and a DataFrame with the parsed data,
               or None if no file is found.
    """
    # Primary pattern for mads_t<num> or rate<num>
    primary_pattern = re.compile(rf"(mads_t|rate){experiment_number:03d}.*\.xls", re.IGNORECASE)
    # Fallback pattern for #<num>
    fallback_pattern = re.compile(rf"#\s*{experiment_number}\.xls", re.IGNORECASE)

    # Search recursively for files
    for file in Path(directory).rglob("*.xls"):
        # Check for primary pattern first
        if primary_pattern.search(file.name):
            try:
                return parse_file(file, sheet_name)
            except Exception as e:
                print(f"Error parsing file {file}: {e}")
                return None

    # If no file is found using the primary pattern, use the fallback pattern
    for file in Path(directory).rglob("*.xls"):
        if fallback_pattern.search(file.name):
            try:
                return parse_file(file, sheet_name)
            except Exception as e:
                print(f"Error parsing file {file}: {e}")
                return None

    # If no file is found
    print(f"No file found for experiment number {experiment_number}.")
    return None


def parse_file(file, sheet_name=None):
    """
    Parses the specified .xls file.

    Args:
        file (Path): Path to the file.
        sheet_name (str, optional): Name of the sheet to parse. If None, the first sheet is used.

    Returns:
        tuple: A tuple containing the filename and a DataFrame with the parsed data.
    """
    # Read the .xls file without headers
    if sheet_name:
        # Read the specified sheet
        data = pd.read_excel(file, sheet_name=sheet_name, header=None)
    else:
        # Read the entire file and default to the first sheet
        df = pd.read_excel(file, sheet_name=None, header=None)  # Load all sheets
        print(f"Available sheets: {list(df.keys())}")
        first_sheet_name = list(df.keys())[0]
        data = df[first_sheet_name]

    return file.name, data


def find_header_row(data, search_strings=None):
    """
    Locates the header row in a given DataFrame by searching for any of the specified strings
    with exact matches, disregarding capitalization.

    Args:
        data (pd.DataFrame): The DataFrame to search through.
        search_strings (list of str): List of strings to search for in the cells.

    Returns:
        int: Index of the header row, or None if no match is found.
    """
    if search_strings is None:
        search_strings = ["[Enz]", "Kuv.", "[Enz] mmol/l", "sub"]  # Default strings to search for

    # Convert the list of search strings to lowercase for case-insensitive matching
    search_strings = [s.lower() for s in search_strings]

    # Iterate through rows of the DataFrame
    for i, row in data.iterrows():
        # Check if any cell matches any search string exactly (case-insensitive)
        if any(str(cell).strip().lower() in search_strings for cell in row):
            return i

    # If no match is found
    print(f"Header row containing one of {search_strings} not found.")
    return None


def find_numeric_values_below_header(data, header_row, sample_num):
    """
    Searches the header row for specific strings (case-insensitive), and for each match,
    finds the first 'sample_num' numeric values below the header row in the corresponding columns.
    Special fallback for [buf] searches the row below the header for strings containing 'buf'.
    Maps the found values to standardized keys: [enz], [h2o2], [sub], and [buf].
    If multiple matches exist, selects values from the column with the largest first numeric value.
    The values are rounded to 3 decimal places. For [buf], adjusts values using (number_in_cell / 2) * 100.

    Args:
        data (pd.DataFrame): The DataFrame containing the data.
        header_row (int): The index of the header row.
        sample_num (int): The number of numeric values to include.

    Returns:
        dict: A dictionary with standardized keys ([enz], [h2o2], [sub], [buf]) and lists of rounded
              numeric values as values.
    """
    results = {"[enz]": None, "[h2o2]": None, "[sub]": None, "[buf]": None}

    # Mapping of standardized keys to search terms
    search_mapping = {
        "[enz]": ["[enz]", "enz"],
        "[h2o2]": ["[h2o2]", "h2o2"],
        "[sub]": ["[sub]", "sub"],
        "[buf]": ["[buf]", "[buffer]"]
    }

    # Extract the header row
    header = data.iloc[header_row]
    row_below_header = data.iloc[header_row + 1]

    # Iterate over the columns in the header
    for col in data.columns:
        cell_value = str(header[col]).lower()  # Convert to lowercase for case-insensitive comparison

        # Check if the cell contains any of the search terms
        for standardized_key, search_terms in search_mapping.items():
            if any(term in cell_value for term in search_terms):
                # Find numeric values below the header
                numeric_values = []
                for value in data[col].iloc[header_row + 1:]:
                    if isinstance(value, (int, float)):
                        numeric_values.append(round(value, 3))
                    if len(numeric_values) == sample_num:
                        break

                if numeric_values:
                    # Store the values, prioritizing the column with the largest first value
                    if results[standardized_key] is None or numeric_values[0] > results[standardized_key][0]:
                        results[standardized_key] = numeric_values

    # Special fallback for [buf]
    if results["[buf]"] is None:
        for col in data.columns:
            cell_value_below = str(row_below_header[col]).lower()  # Row below header
            if "buf" in cell_value_below:
                numeric_values = []
                for value in data[col].iloc[header_row + 2:]:  # Search for numeric values below
                    if isinstance(value, (int, float)):
                        adjusted_value = round((value / 2) * 100, 3)  # Here we assume the buffer is 100 mM
                        numeric_values.append(adjusted_value)
                    if len(numeric_values) == sample_num:
                        break
                if numeric_values:
                    results["[buf]"] = numeric_values
                    break

    # Set any None results to a list of zeros
    for key in results:
        if results[key] is None:
            results[key] = [0] * sample_num

    return results


def _parse_pH_cell(value):
    """
    Reads a pH out of a cell that may be a number or a written range.

    Exp 131 records its buffer pH as the string "8.88-9.07" -- the drift across
    the run rather than a single reading -- so a range is taken at its midpoint.

    Args:
        value: The cell contents.

    Returns:
        float: The pH, or None if the cell holds neither a number nor a range.
    """
    if isinstance(value, (int, float)):
        return None if (isinstance(value, float) and np.isnan(value)) else round(float(value), 2)
    if isinstance(value, str):
        numbers = re.findall(r"\d+(?:[.,]\d+)?", value)
        parsed = [float(n.replace(",", ".")) for n in numbers]
        parsed = [n for n in parsed if 0 < n < 14]
        if len(parsed) == 2:
            return round(sum(parsed) / 2, 2)
        if len(parsed) == 1:
            return round(parsed[0], 2)
    return None


def find_pH_value_in_range(data, row_range, col_range, filename=None):
    """
    Finds the pH value beside a cell labelling it, within a range of rows and
    columns, falling back to the filename patterns "pH=XX.XX" or "pH_XX,XX".

    A sheet may carry more than one pH. Exps 127-131 hold two buffer blocks --
    the pyrophosphate one actually used, labelled "buffer pH", and an unused
    phosphate block whose bare "pH" cell reads 7.29 in all five sheets. Matching
    only the bare label gave all five the same wrong pH and flattened a real
    series spanning 6.94 to 9.0, an error of up to 67x in [HOO-]. So the more
    specific label is searched first. See DATA_VERIFICATION.md 2026-08-31.

    Args:
        data (pd.DataFrame): The DataFrame containing the data.
        row_range (tuple): A tuple specifying the start and end row indices (inclusive).
        col_range (tuple): A tuple specifying the start and end column indices (inclusive).
        filename (str, optional): The filename to parse the pH from as a fallback.

    Returns:
        float: The pH value if found and numeric, or None if not found or is NaN.
    """
    # Extract the subrange of the DataFrame to search
    start_row, end_row = row_range
    start_col, end_col = col_range
    search_area = data.iloc[start_row:end_row + 1, start_col:end_col + 1]

    # Most specific label first: a bare "pH" may belong to a block that was
    # never used, while "buffer pH" names the solution that went in the cuvette.
    for pattern in (r"buffer\s*pH", r"\bpH\b"):
        for row_index in range(search_area.shape[0]):
            for col_index in range(search_area.shape[1] - 1):
                cell = search_area.iloc[row_index, col_index]
                if not (isinstance(cell, str)
                        and re.fullmatch(pattern, cell.strip(), re.IGNORECASE)):
                    continue
                parsed = _parse_pH_cell(search_area.iloc[row_index, col_index + 1])
                if parsed is not None:
                    return parsed

    for row_index in range(search_area.shape[0]):
        for col_index in range(search_area.shape[1] - 1):  # Stop before the last column
            cell = search_area.iloc[row_index, col_index]
            # Match "pH" as a whole word
            if isinstance(cell, str) and re.fullmatch(r"\bpH\b", cell, re.IGNORECASE):
                # Check the cell to the right
                next_cell = search_area.iloc[row_index, col_index + 1]

                # Handle numeric types directly
                if isinstance(next_cell, float) and np.isnan(next_cell):
                    print("Adjacent cell contains NaN.")
                    return None
                if isinstance(next_cell, (int, float)):
                    return round(next_cell, 2)

                print(f"Adjacent value is not numeric: {next_cell}")
                print(filename)
                return None

    # If 'pH' is not found, attempt to parse from the filename
    if filename:
        # First pattern: pH=XX.XX
        match = re.search(r"pH=([\d.]+)", filename, re.IGNORECASE)
        if match:
            return round(float(match.group(1)), 2)

        # Second pattern: pH_XX,XX (comma as decimal separator)
        match = re.search(r"pH_([\d,]+)", filename, re.IGNORECASE)
        if match:
            pH_value = match.group(1).replace(",", ".")  # Convert comma to dot for float conversion
            return round(float(pH_value), 3)

        print("pH value not found in filename.")

    # If neither DataFrame nor filename contains pH information
    print("pH string not found in the specified range or filename.")
    return None


def find_temperature_value_in_range(data, row_range, col_range, filename=None):
    """
    Finds the temperature value located in a cell to the immediate right of a cell containing the string "T",
    within a specified range of rows and columns. If not found, searches for any cell in the range containing
    a string like "XX C" and extracts the numeric temperature value. If still not found, attempts to parse the
    temperature from the filename using the pattern "t=XX".

    Args:
        data (pd.DataFrame): The DataFrame containing the data.
        row_range (tuple): A tuple specifying the start and end row indices (inclusive).
        col_range (tuple): A tuple specifying the start and end column indices (inclusive).
        filename (str, optional): The filename to parse the temperature from as a fallback.

    Returns:
        float: The temperature value if found and numeric, or None if not found.
    """
    # Extract the subrange of the DataFrame to search
    start_row, end_row = row_range
    start_col, end_col = col_range
    search_area = data.iloc[start_row:end_row + 1, start_col:end_col + 1]

    # First attempt: Find "T" and get the adjacent cell
    for row_index in range(search_area.shape[0]):
        for col_index in range(search_area.shape[1] - 1):  # Stop before the last column
            cell = search_area.iloc[row_index, col_index]
            if isinstance(cell, str) and re.fullmatch(r"\bT\b", cell, re.IGNORECASE):
                # Check the cell to the right
                next_cell = search_area.iloc[row_index, col_index + 1]

                # Handle numeric types directly
                if isinstance(next_cell, (int, float)):
                    return round(next_cell, 3)

                print(f"Adjacent value is not numeric: {next_cell}")
                return None

    # Fallback: Search for cells with "XX C" pattern directly in the range
    for row_index in range(search_area.shape[0]):
        for col_index in range(search_area.shape[1]):
            cell = search_area.iloc[row_index, col_index]
            if isinstance(cell, str):
                match = re.search(r"(\d+(\.\d+)?)\s*C", cell, re.IGNORECASE)
                if match:
                    return int(match.group(1))

    # Fallback: Attempt to parse from the filename
    if filename:
        match = re.search(r"t=(\d+)", filename, re.IGNORECASE)
        if match:
            return int(match.group(1))

    # If no temperature information is found
    print("Temperature value not found in the specified range or filename.")
    return None


def find_buffer_type(data, row_range, col_range, filename=None, buffer_map=None):
    """
    Finds the buffer type within a specified range of rows and columns in a DataFrame.
    If not found in the DataFrame, attempts to parse the buffer type from the filename.
    Maps found buffer types to standardized names using a buffer_map.

    Args:
        data (pd.DataFrame): The DataFrame containing the data.
        row_range (tuple): A tuple specifying the start and end row indices (inclusive).
        col_range (tuple): A tuple specifying the start and end column indices (inclusive).
        filename (str, optional): The filename to parse the buffer type from as a fallback.
        buffer_map (dict, optional): Mapping of buffer type keywords to standardized names.

    Returns:
        str: The standardized buffer type if found, or None if not found.
    """
    if buffer_map is None:
        buffer_map = {
            "boric": "Boric",
            "pyrophosphate": "Pyrophosphate",
            "phosphat": "Phosphate",
            "phosphate": "Phosphate",
            "Na4P2O7*10H2O": "Pyrophosphate",
            "NaH2PO4*2H2O": "Phosphate",   # monosodium phosphate, NOT pyrophosphate
            "carbonate": "Carbonate",
            "CO3": "Carbonate"
        }

    # Create a set of all keywords to search for (case-insensitive)
    keywords = {key.lower(): value for key, value in buffer_map.items()}

    # Extract the subrange of the DataFrame to search
    start_row, end_row = row_range
    start_col, end_col = col_range
    search_area = data.iloc[start_row:end_row + 1, start_col:end_col + 1]

    # Search the DataFrame for buffer types
    for row_index in range(search_area.shape[0]):
        for col_index in range(search_area.shape[1]):
            cell = search_area.iloc[row_index, col_index]
            if isinstance(cell, str):
                cell_lower = cell.lower()
                for keyword, standardized_name in keywords.items():
                    if keyword in cell_lower:
                        return standardized_name

    # If no buffer type is found, attempt to parse from the filename
    if filename:
        filename_lower = filename.lower()
        for keyword, standardized_name in keywords.items():
            if keyword in filename_lower:
                return standardized_name

    # If neither DataFrame nor filename contains buffer type information
    print("Buffer type not found in the specified range or filename.")
    print(f"Filename: {filename}")
    return None


def find_substrate_type(data, row_range, col_range, filename=None):
    """
    Finds the substrate type within a specified range of rows and columns in a DataFrame.
    If not found in the DataFrame, attempts to parse the substrate type from the filename.
    Maps recognized substrate types to properly formatted values.

    Args:
        data (pd.DataFrame): The DataFrame containing the data.
        row_range (tuple): A tuple specifying the start and end row indices (inclusive).
        col_range (tuple): A tuple specifying the start and end column indices (inclusive).
        filename (str, optional): The filename to parse the substrate type from as a fallback.

    Returns:
        str: The properly formatted substrate type if found, or an empty string if not found.
    """
    # Substrate types and their standardized mappings
    substrate_mapping = {
        "bnoh": "BnOH",
        "benzylalkohol": "BnOH",
        "4ome-bnoh": "4OMe-BnOH",
        "4-meo-bnoh": "4OMe-BnOH",
        "4-methoxy-benzylalkohol": "4OMe-BnOH"
    }

    # Convert mapping keys to lowercase for case-insensitive comparison
    search_terms = {key.lower(): value for key, value in substrate_mapping.items()}

    # Extract the subrange of the DataFrame to search
    start_row, end_row = row_range
    start_col, end_col = col_range
    search_area = data.iloc[start_row:end_row + 1, start_col:end_col + 1]

    # Search the DataFrame for substrate types
    for row_index in range(search_area.shape[0]):
        for col_index in range(search_area.shape[1]):
            cell = search_area.iloc[row_index, col_index]
            if isinstance(cell, str):
                cell_lower = cell.lower()
                for term, standardized in search_terms.items():
                    # Match exact, start-of-string followed by whitespace, or whitespace-enclosed matches
                    if re.search(rf"(^|\s){re.escape(term)}(\s|$)", cell_lower):
                        return standardized

    # If no substrate type is found, attempt to parse from the filename
    if filename:
        filename_lower = filename.lower()
        for term, standardized in search_terms.items():
            if re.search(rf"(^|\s){re.escape(term)}(\s|$)", filename_lower):
                return standardized

    # If neither DataFrame nor filename contains substrate type information
    print("Substrate type not found in the specified range or filename.")
    print(f"Filename: {filename}")
    return ""


def populate_experimental_data_from_directory(directory, sheet_name='Sheet1'):
    """
    Populates a DataFrame with experimental data extracted from all experiment files in the directory.
    Additionally returns a list of all parsed_data objects.

    Args:
        directory (str): Directory containing both the experiment data files (.txt) and associated .xls files.
        sheet_name (str): Name of the sheet to parse in the .xls files.

    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: DataFrame containing the experimental data.
            - list: List of parsed_data objects from all text files.
    """
    experiment_rows = []
    parsed_data_list = []

    text_files = [
        os.path.join(directory, file)
        for file in os.listdir(directory)
        if file.endswith(".txt")
    ]

    for file_path in text_files:
        try:
            parsed_data = parse_experiment_data(file_path)
            parsed_data_list.append(parsed_data)
            experiment_number = parsed_data["num"]

            file_name, experiment_data = find_and_parse_experiment_file(
                experiment_number=experiment_number,
                directory=directory,
                sheet_name=sheet_name
            )

            header_row = find_header_row(experiment_data)
            initial_cons = find_numeric_values_below_header(
                experiment_data, header_row, sample_num=len(parsed_data['samples'])
            )
            pH_value = find_pH_value_in_range(experiment_data, (0, 100), (0, 100), filename=file_name)
            temperature = find_temperature_value_in_range(experiment_data, (0, 100), (0, 100), filename=file_name)
            buffer_type = find_buffer_type(experiment_data, (0, 100), (0, 100), filename=file_name)
            substrate_type = find_substrate_type(experiment_data, (0, 100), (0, 100), filename=file_name)

            abs_value = SUBSTRATE_PROPERTIES.get(substrate_type, {}).get("abs", None)
            e_value = SUBSTRATE_PROPERTIES.get(substrate_type, {}).get("e", None)

            for i, sample_name in enumerate(parsed_data['samples']):
                cons = apply_experiment_corrections(experiment_number, i, {
                    "[enz]": initial_cons["[enz]"][i] if initial_cons["[enz]"] else 0,
                    "[buf]": initial_cons["[buf]"][i] if initial_cons["[buf]"] else 0,
                    "[h2o2]": initial_cons["[h2o2]"][i] if initial_cons["[h2o2]"] else 0,
                    "[sub]": initial_cons["[sub]"][i] if initial_cons["[sub]"] else 0,
                })
                row = {
                    "experiment": experiment_number,
                    "sample": i + 1,
                    "substrate": substrate_type,
                    "abs": abs_value,
                    "e": e_value,
                    "buffer": buffer_type,
                    "pH": pH_value,
                    "T": temperature,
                    **cons,
                }
                experiment_rows.append(row)

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            continue

    experiment_df = pd.DataFrame(experiment_rows)

    return experiment_df, parsed_data_list


def load_experiment(experiment_number, directory="data/data", sheet_name="Sheet1"):
    """
    Loads one experiment's metadata and raw time series together, for tools
    (e.g. plot_kinetics.py) that want a single experiment rather than the
    whole-dataset DataFrame. Uses exactly the same extraction calls as
    populate_experimental_data_from_directory, just scoped to one experiment.

    Returns:
        dict: {
            "experiment": int,
            "xls_file": str,
            "txt_file": str,
            "substrate": str, "buffer": str, "pH": float, "T": float,
            "samples": [
                {"sample": 1, "sample_name": "Sample001",
                 "[enz]": .., "[buf]": .., "[h2o2]": .., "[sub]": ..,
                 "abs": .., "e": ..,
                 "time": [...], "values": [...]},
                ...
            ],
        }
        or None if the experiment or its files can't be found.
    """
    txt_path = None
    for file in os.listdir(directory):
        if file.endswith(".txt"):
            candidate = os.path.join(directory, file)
            parsed = parse_experiment_data(candidate)
            if parsed is not None and parsed.get("num") == experiment_number:
                txt_path = candidate
                parsed_data = parsed
                break
    if txt_path is None:
        print(f"No .txt file found for experiment {experiment_number}.")
        return None

    file_name, experiment_data = find_and_parse_experiment_file(
        experiment_number=experiment_number, directory=directory, sheet_name=sheet_name
    )
    if experiment_data is None:
        print(f"No .xls file found for experiment {experiment_number}.")
        return None

    header_row = find_header_row(experiment_data)
    n_samples = len(parsed_data["samples"])
    initial_cons = find_numeric_values_below_header(experiment_data, header_row, sample_num=n_samples)
    pH_value = find_pH_value_in_range(experiment_data, (0, 100), (0, 100), filename=file_name)
    temperature = find_temperature_value_in_range(experiment_data, (0, 100), (0, 100), filename=file_name)
    buffer_type = find_buffer_type(experiment_data, (0, 100), (0, 100), filename=file_name)
    substrate_type = find_substrate_type(experiment_data, (0, 100), (0, 100), filename=file_name)
    abs_value = SUBSTRATE_PROPERTIES.get(substrate_type, {}).get("abs")
    e_value = SUBSTRATE_PROPERTIES.get(substrate_type, {}).get("e")

    samples = []
    for i, sample_name in enumerate(parsed_data["samples"]):
        s = parsed_data["samples"][sample_name]
        cons = apply_experiment_corrections(experiment_number, i, {
            "[enz]": initial_cons["[enz]"][i] if initial_cons["[enz]"] else 0,
            "[buf]": initial_cons["[buf]"][i] if initial_cons["[buf]"] else 0,
            "[h2o2]": initial_cons["[h2o2]"][i] if initial_cons["[h2o2]"] else 0,
            "[sub]": initial_cons["[sub]"][i] if initial_cons["[sub]"] else 0,
        })
        samples.append({
            "sample": i + 1,
            "sample_name": sample_name,
            **cons,
            "abs": abs_value,
            "e": e_value,
            "time": s["time"],
            "values": s["values"],
        })

    return {
        "experiment": experiment_number,
        "xls_file": file_name,
        "txt_file": os.path.basename(txt_path),
        "substrate": substrate_type,
        "buffer": buffer_type,
        "pH": pH_value,
        "T": temperature,
        "samples": samples,
    }
