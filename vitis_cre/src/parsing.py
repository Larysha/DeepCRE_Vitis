from __future__ import annotations
import argparse
from enum import Enum, auto
import json
import os
from typing import Any, Dict, Iterable, List, Tuple
from copy import deepcopy
import pandas as pd

from vitis_cre.src.utils import make_absolute_path, read_feature_from_input_dict



class RunInfo:
    """holds information about a single run of a deepCRE script.

    Information is split over general information that is used for all species, and species specific information.
    """
    general_info: Dict[str, Any]
    species_info: List[Dict[str, Any]]
    possible_general_parameters: Dict
    possible_species_parameters: Dict

    def __init__(self, possible_general_parameters: Dict, possible_species_parameters: Dict) -> None:
        """initialize the RunInfo object.

        Args:
            possible_general_parameters (Dict): defines the possible parameters that can be used in the general info dict and their default values.
            possible_species_parameters (Dict): defines the possible parameters that can be used in the species info dict and their default values.
        """
        self.possible_general_parameters = possible_general_parameters
        self.possible_species_parameters = possible_species_parameters

    def set_up_defaults(self) -> Dict[str, Any]:
        """set up the default values for the species info dict based on the general info dict and the possible parameters.

        Returns:
            Dict[str, Any]: dictionary containing the default values for the species info dict.
        """
        #get general defaults
        defaults = self.possible_species_parameters.copy()
        #get values from general columns that can be used for species infos
        applicable_general_inputs = list(set(self.possible_general_parameters.keys()).intersection(set(self.possible_species_parameters.keys())))
        applicable_general_inputs = {key: self.general_info[key] for key in applicable_general_inputs if key in self.general_info.keys()}
        #overwrite the defaults with info from general_info, where applicable
        defaults.update(applicable_general_inputs)
        applicable_general_inputs = defaults
        #remove key/value pairs that dont have a meaningful default or a value from the general_info, so they dont interfere later
        to_delete = [key for key, value in applicable_general_inputs.items() if value is None]
        for key in to_delete:
            del applicable_general_inputs[key]
        return applicable_general_inputs

    def check_general_keys(self, general_dict: Dict[str, Any]) -> None:
        """check if all necessary parameters are present in the general info dict.

        Args:
            general_dict (Dict[str, Any]): dictionary containing the general info.

        Raises:
            ValueError: if any necessary parameter is missing.
        """
        necessary_parameters = list(set(self.possible_general_parameters).difference(set(self.possible_species_parameters)))
        missing = [parameter for parameter in necessary_parameters if parameter not in general_dict.keys()]
        if missing:
            raise ValueError(f"Error parsing general info dict! Input parameters {missing} missing!")

    def check_species_keys(self, specie_dict: Dict[str, Any], missing_species_name_ok: bool = False):
        """check if all necessary parameters are present in the species info dict.

        Args:
            specie_dict (Dict[str, Any]): dictionary containing the species info.
            missing_species_name_ok (bool, optional): if True, the function will not raise an error if only the species_name is missing. Defaults to False.

        Raises:
            ValueError: if any necessary parameter is missing.
        """
        missing = [parameter for parameter in self.possible_species_parameters.keys() if parameter not in specie_dict.keys()]
        if missing_species_name_ok and missing == ["species_name"]:
            return
        if missing:
            raise ValueError(f"Error parsing specie dict! Input parameters {missing} missing!")
    
    def load_chromosomes(self) -> None:
        """load the chromosomes from the csv file if a path is given instead of a list of chromosomes.

        Raises:
            ValueError: if the input for the chromosomes is not a list or a path to a csv file.
        """
        for specie_info in self.species_info:
            if not "chromosomes" in specie_info.keys():
                continue
            if specie_info["chromosomes"] == "":
                specie_info["chromosomes"] = []
                continue
            chromosomes = specie_info["chromosomes"]
            if isinstance(chromosomes, list):
                continue
            elif not isinstance(chromosomes, str):
                raise ValueError(f"input for chomosomes needs to be a path to a csv containing the chromosomes in a single column (type: str), or a list containing the names of the chromosomes (tpye: list of strings). Found type: {type(chromosomes)}.")
            path = chromosomes if os.path.isfile(chromosomes) else make_absolute_path("genome", chromosomes, start_file=__file__)
            chromosomes = pd.read_csv(path, header=None).values.ravel().tolist()
            chromosomes = [str(chrom) for chrom in chromosomes]
            specie_info["chromosomes"] = chromosomes


    def parse_species(self, run_dict: Dict[str, Any]) -> None:
        """parse the species info from the input dictionary.

        Args:
            run_dict (Dict[str, Any]): dictionary that was read from the input file
        """
        #list of parameters in both general and species parameters
        applicable_general_inputs = self.set_up_defaults()
        species_info = []
        for species_dict in run_dict.get("species_data", []):
            curr_specie_dict = {key: read_feature_from_input_dict(species_dict, key) for key in self.possible_species_parameters.keys() if key in species_dict.keys()}
            #get defaults, then overwrite them with the data that was read in for all keys where data was read in
            curr_general_info = applicable_general_inputs.copy()
            curr_general_info.update(curr_specie_dict)
            #make sure all necessary parameters are filled
            self.check_species_keys(curr_general_info)
            species_info.append(curr_general_info)
        if not species_info:
            curr_general_info = applicable_general_inputs.copy()
            self.check_species_keys(curr_general_info, missing_species_name_ok=True)
            species_info.append(curr_general_info)
        self.species_info = species_info
        self.load_chromosomes()
    
    def load_prediction_models(self) -> None:
        """load the prediction models from the csv file if a path is given instead of a list of model names.

        Raises:
            ValueError: if the input for the prediction models is not a list or a path to a csv file.
        """
        if not "prediction_models" in self.general_info.keys():
            return
        prediction_models = self.general_info["prediction_models"]
        if isinstance(prediction_models, list):
            return
        elif not isinstance(prediction_models, str):
            raise ValueError(f"input for prediction models needs to be a path to a csv containing the models in a single column (type: str), or a list containing the names of the chromosomes (type: list of strings). Found type: {type(prediction_models)}.")

        path = prediction_models if os.path.isfile(prediction_models) else make_absolute_path("saved_models", prediction_models, start_file=__file__)
        prediction_models = pd.read_csv(path, header=None).values.ravel().tolist()
        prediction_models = [str(chrom) for chrom in prediction_models]
        self.general_info["prediction_models"] = prediction_models

    def parse_general_inputs(self, run_dict: Dict[str, str]) -> None:
        """parse the general info from the input dictionary.

        Args:
            run_dict (Dict[str, str]): dictionary that was read from the input file
        """
        #load defaults
        defaults = self.possible_general_parameters.copy()
        read_data = {key: read_feature_from_input_dict(run_dict, key) for key in self.possible_general_parameters.keys() if key in run_dict.keys()}
        #overwrite defaults with read data
        defaults.update(read_data)
        general_info = defaults
        #remove empty defaults
        to_delete = [key for key, value in general_info.items() if value is None]
        for key in to_delete:
            del general_info[key]
        if "model_case" in general_info.keys():
            general_info["model_case"] = ModelCase.parse(general_info["model_case"])
        #make check
        self.check_general_keys(general_dict=general_info)
        self.general_info = general_info
        self.load_prediction_models()

    @staticmethod
    def parse(run_dict: Dict[str, Any], possible_general_parameters: Dict, possible_species_parameters: Dict, multiple_species_required_msr: bool = False) -> RunInfo:
        """parse the input dictionary into a RunInfo object.

        Args:
            run_dict (Dict[str, Any]): dictionary that was read from the input file
            possible_general_parameters (Dict): Dict containing the possible parameters for the general info dict and their default values.
            possible_species_parameters (Dict): Dict containing the possible parameters for the species info dict and their default values.
            multiple_species_required_msr (bool, optional): if True, the function will raise an error if less than 2 species are present in the species info dict for msr runs. Defaults to False.

        Raises:
            ValueError: if less than 2 species are present in the species info dict for msr runs.

        Returns:
            RunInfo: RunInfo object containing the parsed information.
        """
        run_info_object = RunInfo(possible_general_parameters=possible_general_parameters, possible_species_parameters=possible_species_parameters)
        run_info_object.parse_general_inputs(run_dict=run_dict)
        # load general info first, and use it as defaults for specific species
        run_info_object.parse_species(run_dict)
        if multiple_species_required_msr and "model_case" in run_info_object.general_info.keys() and run_info_object.general_info["model_case"] == ModelCase.MSR and len(run_info_object.species_info) < 2:
            raise ValueError(f"Need at least 2 species for MSR training! Only found {len(run_info_object.species_info)}.")
        for species_key in possible_species_parameters.keys():
            if species_key in run_info_object.general_info.keys():
                del run_info_object.general_info[species_key]
        return run_info_object
    
    def is_msr(self) -> bool:
        """check if the run is a MSR run.

        Returns:
            bool: True if the run is a MSR run, False otherwise.
        """
        if "model_case" in self.general_info.keys():
            return self.general_info["model_case"] == ModelCase.MSR
        else:
            return False
    
    def __str__(self) -> str:
        """return a string representation of the RunInfo object.

        Returns:
            str: string representation of the RunInfo object based on the json format.
        """
        if "model_case" in self.general_info.keys():
            self.general_info["model_case"] = str(self.general_info["model_case"])
        specie_info_json = json.dumps(self.species_info, indent=2)
        gen_info_json = json.dumps(self.general_info, indent=2)
        if "model_case" in self.general_info.keys():
            self.general_info["model_case"] = ModelCase.parse(self.general_info["model_case"])
        result = "RunInfo(\n  General Info{"
        for gen_info in gen_info_json.split("\n"):
            if gen_info in ["{", "}"]:
                continue
            result += "\n  " + gen_info
        result += "\n  },\n  Species info["
        for species_info in specie_info_json.split("\n"):
            if species_info in ["[", "]"]:
                continue
            result += "\n  " + species_info
        result += "\n  ]\n)"
        return result


class ParsedInputs:
    """holds information about multiple runs of a deepCRE script.

    Holds possible parameters for general and species info, as well as their defaults, since they are consistent for scripts.
    """
    run_infos: List[RunInfo]
    possible_general_parameters: Dict
    possible_species_parameters: Dict

    def __init__(self, possible_general_parameters: Dict, possible_species_parameters: Dict):
        """initialize the ParsedInputs object.

        Args:
            possible_general_parameters (Dict): possible parameters for the general info dict and their default values.
            possible_species_parameters (Dict): possible parameters for the species info dict and their default values.
        """
        self.run_infos = []
        self.possible_general_parameters = possible_general_parameters
        self.possible_species_parameters = possible_species_parameters

    @staticmethod
    def parse(json_file_name: str, possible_general_parameters: Dict, possible_species_parameters: Dict, allow_multiple_species: bool = True, multiple_species_required_msr: bool = False) -> Tuple[ParsedInputs, List[Tuple[str, int, Exception]], int]:
        """parse the input json file into a ParsedInputs object.

        Args:
            json_file_name (str): path to the json file containing the input data.
            possible_general_parameters (Dict): possible_general_parameters (Dict): _description_
            possible_species_parameters (Dict): possible_species_parameters (Dict): _description_
            allow_multiple_species (bool, optional): Determines whether it is allowed to have multiple species in a single run. Defaults to True.
            multiple_species_required_msr (bool, optional): Will fail msr runs with less than 2 species. Defaults to False.

        Returns:
            Tuple[ParsedInputs, List[Tuple[str, int, Exception]], int]: Tuple containing the ParsedInputs object, a list of failed parsings, and the total number of runs.
        """
        json_file_name = json_file_name if os.path.isfile(json_file_name) else make_absolute_path(json_file_name, __file__)
        failed_parsings = []
        parsed_object = ParsedInputs(possible_general_parameters=possible_general_parameters, possible_species_parameters=possible_species_parameters)
        with open(json_file_name, "r") as f:
            input_list = json.load(f)
        for i, run_dict in enumerate(input_list):
            try:
                curr_run_info = RunInfo.parse(run_dict, possible_general_parameters=parsed_object.possible_general_parameters, possible_species_parameters=parsed_object.possible_species_parameters, multiple_species_required_msr=multiple_species_required_msr)
                if len(curr_run_info.species_info) > 1 and not allow_multiple_species:
                    raise ValueError(f"Only one species per run allowed for running this script! Found {len(curr_run_info.species_info)}!")
                parsed_object.run_infos.append(curr_run_info)
            except Exception as e:
                print(f"error reading input run number {i + 1}.")
                print(f"error message is: \"{e}\"")
                print(f"the dictionary that was loaded for the run is the following:")
                print(f"{json.dumps(run_dict, indent=2)}")
                failed_parsings.append(("error during parsing!", i, e))
        return parsed_object, failed_parsings, len(input_list)
    
    def replace_both(self) -> ParsedInputs:
        """replace runs with model_case "both" with two runs, one for SSR and one for SSC.

        Returns:
            ParsedInputs: new ParsedInputs object with the replaced runs.
        """
        new_object = ParsedInputs(possible_species_parameters=self.possible_species_parameters, possible_general_parameters=self.possible_general_parameters)
        for info in self.run_infos:
            if "model_case" in info.general_info.keys() and info.general_info["model_case"] == ModelCase.BOTH:
                info.general_info["model_case"] = ModelCase.SSR
                ssc_version = deepcopy(info)
                ssc_version.general_info["model_case"] = ModelCase.SSC
                new_object.run_infos.append(info)
                new_object.run_infos.append(ssc_version)
            else:
                new_object.run_infos.append(info)
        return new_object
    
    def __iter__(self) -> Iterable[RunInfo]:
        """return an iterable over the RunInfo objects.

        Returns:
            Iterable[RunInfo]: iterable over the RunInfo objects.
        """
        return (info for info in self.run_infos)
    
    def __str__(self) -> str:
        """return a string representation of the ParsedInputs object.

        Returns:
            str: string representation of the ParsedInputs object based on the json format.
        """
        result = "ParsedTrainingInputs["
        for info in self.run_infos:
            for info_line in str(info).split("\n"):
                result += "\n  " + info_line
        result += "\n]"
        return result
    
    def __len__(self) -> int:
        """return the number of runs in the ParsedInputs object.

        Returns:
            int: number of runs in the ParsedInputs object.
        """
        return len(self.run_infos)
            


class ModelCase(Enum):
    """Enum to represent the different model cases that can be used in a deepCRE run.
    """
    MSR = auto()
    SSR = auto()
    SSC = auto()
    BOTH = auto()

    @staticmethod
    def parse(input_string: str) -> ModelCase:
        """parse a string into a ModelCase object.

        Args:
            input_string (str): string to be parsed.

        Raises:
            ValueError: if the input string is not recognized.

        Returns:
            ModelCase: ModelCase object representing the input string.
        """
        if input_string.lower() == "msr":
            return ModelCase.MSR
        elif input_string.lower() == "ssr":
            return ModelCase.SSR
        elif input_string.lower() == "ssc":
            return ModelCase.SSC
        elif input_string.lower() == "both":
            return ModelCase.BOTH
        else:
            raise ValueError(f"model case \"{input_string}\" not recognized!")
    
    def __str__(self) -> str:
        """return a string representation of the ModelCase object.

        Raises:
            NotImplementedError: if the string method was not implemented for this Variant yet!

        Returns:
            str: string representation of the ModelCase object.
        """
        if self == ModelCase.MSR:
            return "msr"
        elif self == ModelCase.SSR:
            return "ssr"
        elif self == ModelCase.SSC:
            return "ssc"
        elif self == ModelCase.BOTH:
            return "both"
        else:
            raise NotImplementedError("string method was not implemented for this Variant yet!")
    

def identity(obj: Any) -> Any:
    """identity function that returns the input object.

    Args:
        obj (Any): input object.

    Returns:
        Any: input object.
    """
    return obj


def convert_named_csv(csv_path: str, output_json_path: str, **command_line_args) -> None:
    """converts a csv file with named columns to a json file.

    Args:
        csv_path (str): path to the csv file to be converted
        output_json_path (str): path to the save location of the converted json output
        command_line_args (_type_): additional command line arguments that were used with the script.
    """
    df = pd.read_csv(csv_path, header=0)
    # float and int from pandas are not json serializable, so we need to convert them to python types
    # step 1: get the types of the columns and save which conversion function to use
    types = df.dtypes
    conversions = {}
    for i, (col, dtype) in enumerate(types.items()):
        if str(dtype).startswith("int"):
            conversions[col] = int
        elif str(dtype).startswith("float"):
            conversions[col] = float
        else:
            conversions[col] = identity
    json_compatible_inputs = []
    for i, row in df.iterrows():
        copied = command_line_args.copy()
        # step 2: convert the values in the row to python types
        read_in = {key: conversions[key](row[key]) for key in row.keys()}
        copied.update(read_in)
        json_compatible_inputs.append(copied)
    with open(output_json_path, "w") as f:
        json.dump(json_compatible_inputs, f, indent=2)


def convert_csv_to_json(csv_path: str, output_json_path: str, parsing_mode: str = "training", **command_line_args) -> None:
    """ converts existing csv input files for the older version of the code to the new json format.

    Args:
        csv_path (str): path to the csv file to be converted
        output_json_path (str): path to the save location of the converted json output
        parsing_mode (str, optional): mode in which the csv file'is supposed to be parsed. "training" for old ssr training files, "prediction" for old ssr
            predicition, interpretation or motif extraction files and "named" for csv files with named columns. Defaults to "training".
    """
    if parsing_mode == "named":
        convert_named_csv(csv_path, output_json_path, **command_line_args)
        return
    df = pd.read_csv(csv_path, header=None)
    json_compatible_inputs = []
    for i, (row) in df.iterrows():
        copied = command_line_args.copy()
        if parsing_mode == "training":
            genome, gtf, tpm, output, chroms, p_key = row
            read_in = {
                "genome": genome,
                "annotation": gtf,
                "targets": tpm,
                "output_name": output,
                "chromosomes": chroms,
                "pickle_key": p_key,
            }
            copied.update(read_in)
            json_compatible_inputs.append(copied)
        elif parsing_mode == "prediction":
            genome, gtf, tpm, output, chroms = row
            read_in = {
                "genome": genome,
                "annotation": gtf,
                "targets": tpm,
                "output_name": output,
                "chromosomes": chroms,
            }
            copied.update(read_in)
            json_compatible_inputs.append(copied)
        elif parsing_mode != "named":
            raise ValueError(f"parsing mode \"{parsing_mode}\" not recognized! Must be one of \"training\", \"prediction\" or \"named\".")

    with open(output_json_path, "w") as f:
        json.dump(json_compatible_inputs, f, indent=2)


def parse_args() -> Tuple[str, str, str, Dict[str, str]]:
    """parse command line arguments.

    will also convert the command line arguments after \"--command_line_args\" into a dictionary, interpreting them as pairs of key and value.

    Returns:
        Tuple[str, str, Dict[str, str]]: path to the csv file to be converted, path to the save location of the converted json output, dictionary containing the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_path", help="path to the csv file to be converted", required=True, type=str)
    parser.add_argument("--output_path", help="path to the save location of the converted json output", required=True, type=str)
    parser.add_argument("--parsing_mode", help="mode in which the csv file'is supposed to be parsed. \"training\" for old ssr training files, \"prediction\" for old ssr predicition," +
                        " interpretation or motif extraction files and \"named\" for csv files with named columns.", required=True, type=str)
    parser.add_argument("--command_line_args", help="command line arguments that would have been used with this script", required=True, type=str, nargs="+")
    args = parser.parse_args()
    cmd_line_list = args.command_line_args
    cmd_line_args = {}
    for i in range(0, len(cmd_line_list), 2):
        cmd_line_args[cmd_line_list[i]] = cmd_line_list[i+1]
    for  key, value in cmd_line_args.items():
        try:
            cmd_line_args[key] = int(value)
        except ValueError:
            try:
                cmd_line_args[key] = float(value)
            except ValueError:
                pass
    return args.csv_path, args.output_path, args.parsing_mode, cmd_line_args

if __name__ == "__main__":
    csv, out, mode, args = parse_args()
    convert_csv_to_json(csv_path=csv, output_json_path=out, parsing_mode=mode, **args)

