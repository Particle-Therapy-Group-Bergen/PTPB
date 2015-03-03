#!/usr/bin/env python

# This script provides a driver for the processPatients.m function,
# i.e. provides convenient shell command line tool for processing DVH patient
# data.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    Particle Therapy Project Bergen (PTPB) - tools and models for research in
#    cancer therapy using particle beams.
#
#    Copyright (C) 2015 Particle Therapy Group Bergen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Artur Szostak <artursz@iafrica.com>

from __future__ import print_function
import os
import sys
import argparse
import collections
import subprocess
import textwrap


class RunError(Exception):
    """Exception class for run time errors handled by this script."""
    pass


class SetUnique(argparse.Action):
    """
    A parser action class that ensures an argument was only set once on the
    command line.
    """

    def __init__(self, *args, **kwargs):
        argparse.Action.__init__(self, *args, **kwargs)
        self.alreadyset = False

    def __call__(self, parser, namespace, value, option_string=None):
        if self.alreadyset:
            raise argparse.ArgumentError(self, "can only be given once")
        setattr(namespace, self.dest, value)
        self.alreadyset = True


class UncertaintyModel(object):
    """Base class for uncertainty models."""

    def get_name(self):
        """Return the name of uncertainty model for the Matlab code."""
        raise NotImplementedError

    def get_range(self):
        """Return the uncertainty range compatible with DVH data points."""
        raise NotImplementedError

    def as_matlab(self):
        """Return the uncertainty model encoded as a Matlab code snipet."""
        raise NotImplementedError


class TwoValueRange(UncertaintyModel):
    """Two parameter uncertainty model defining a range [low .. high]."""

    def __init__(self, *args):
        if len(args) == 1:
            self.low = -args[0]
            self.high = args[0]
        elif len(args) == 2:
            if args[0] > args[1]:
                msg = "The second argument '{0}' must be bigger than the" \
                      " first '{1}'.".format(args[1], args[0])
                raise SyntaxError(msg)
            self.low = args[0]
            self.high = args[1]
        else:
            msg = "Maximum of 2 arguments supported for the {0}" \
                  " uncertainty model.".format(self.__class__.__name__)
            raise SyntaxError(msg)

    def as_matlab(self):
        paramstr = "{{" + str(self.low) + ", " + str(self.high) + "}}"
        return "struct('uncertainty_model', '{1}', 'params', {2})".format(
                    self.get_name(), paramstr)


class ThreeValueRange(UncertaintyModel):
    """
    Three parameter uncertainty model defining a range [low .. high] with a
    mode within that range.
    """

    def __init__(self, *args):
        if len(args) == 1:
            self.low = -args[0]
            self.mode = 0
            self.high = args[0]
        elif len(args) == 2:
            if args[0] > args[1]:
                msg = "The second argument '{0}' must be bigger than the" \
                      " first '{1}'.".format(args[1], args[0])
                raise SyntaxError(msg)
            self.low = args[0]
            self.mode = (args[0] + args[1]) * 0.5
            self.high = args[1]
        elif len(args) == 3:
            if args[0] > args[1]:
                msg = "The second argument '{0}' must be bigger than the" \
                      " first '{1}'.".format(args[1], args[0])
                raise SyntaxError(msg)
            if args[1] > args[2]:
                msg = "The third argument '{0}' must be bigger than the" \
                      " second '{1}'.".format(args[2], args[1])
                raise SyntaxError(msg)
            self.low = args[0]
            self.mode = args[1]
            self.high = args[2]
        else:
            msg = "Maximum of 3 arguments supported for the {0}" \
                  " uncertainty model.".format(self.__class__.__name__)
            raise SyntaxError(msg)

    def as_matlab(self):
        paramstr = "{{" + str(self.low) + ", " + str(self.mode) + ", " \
                   + str(self.high) + "}}"
        return "struct('uncertainty_model', '{0}', 'params', {1})".format(
                    self.get_name(), paramstr)


class Delta(UncertaintyModel):
    """Single value delta function distribution."""

    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return "Delta at {0}".format(self.value)

    def get_name(self):
        return "delta"

    def get_range(self):
        return 0

    def as_matlab(self):
        paramstr = "{{" + str(self.value) + "}}"
        return "struct('uncertainty_model', '{0}', 'params', {1})".format(
                    self.get_name(), paramstr)


class DoubleDelta(TwoValueRange):
    """Two value distribution formed by two delta functions."""

    def __repr__(self):
        return "Double delta at {0} and {1}".format(self.low, self.high)

    def get_name(self):
        return "double_delta"


class Box(TwoValueRange):
    """Box distribution."""

    def __repr__(self):
        return "Box uncertainty [{0} .. {1}]".format(self.low, self.high)

    def get_name(self):
        return "box"

    def get_range(self):
        return (self.high - self.low) * 0.5


class Box95(Box):
    """Box distribution defined by a 95% confidence interval range."""

    def __repr__(self):
        return "Box uncertainty with 95% confidence interval" \
               " [{0} .. {1}]".format(self.low, self.high)

    def get_name(self):
        return "box95"

    def get_range(self):
        return (self.high - self.low) * 1.0/0.95 * 0.5


class Triangle(ThreeValueRange):
    """Triangle distribution."""

    def __repr__(self):
        return "Triangle uncertainty [{0} .. {1}] and mode at" \
               " {2}".format(self.low, self.high, self.mode)

    def get_name(self):
        return "triangle"

    def get_range(self):
        return (self.high - self.low) * 0.5


class Triangle95(Triangle):
    """
    Triangle distribution defined by a 95% confidence interval range using the
    mean of the distribution as a parameter.
    """

    def __init__(self, *args):
        super(Triangle, self).__init__(*args)
        self.mean = self.mode
        del self.mode

    def __repr__(self):
        return "Triangle uncertainty with 95% confidence interval [{0} .." \
               " {1}] and mean at {2}".format(self.low, self.high, self.mean)

    def get_name(self):
        return "triangle95"

    def get_range(self):
        return (self.high - self.low) * 1.0/0.95 * 0.5

    def as_matlab(self):
        paramstr = "{{" + str(self.low) + ", " + str(self.mean) + ", " \
                   + str(self.high) + "}}"
        return "struct('uncertainty_model', '{0}', 'params', {1})".format(
                    self.get_name(), paramstr)


class Triangle95mode(Triangle):
    """
    Triangle distribution defined by a 95% confidence interval range using the
    mode of the distribution as a parameter.
    """

    def __repr__(self):
        return "Triangle uncertainty with 95% confidence interval [{0} .." \
               " {1}] and mode at {2}".format(self.low, self.high, self.mode)

    def get_name(self):
        return "triangle95mode"

    def get_range(self):
        return (self.high - self.low) * 1.0/0.95 * 0.5


class Gaus(UncertaintyModel):
    """Gaussian distribution defined by a mean and sigma parameter."""

    def __init__(self, *args):
        if len(args) == 1:
            self.mean = 0
            self.sigma = args[0]
        elif len(args) == 2:
            self.mean = args[0]
            self.sigma = args[1]
        else:
            msg = "Maximum of 2 arguments supported for the {0}" \
                  " uncertainty model.".format(self.__class__.__name__)
            raise SyntaxError(msg)

    def __repr__(self):
        return "Gaussian uncertainty mean = {0}, standard deviation =" \
               " {1}".format(self.mean, self.sigma)

    def get_name(self):
        return "gaus"

    def as_matlab(self):
        paramstr = "{{" + str(self.mean) + ", " + str(self.sigma) + "}}"
        return "struct('uncertainty_model', '{0}', 'params', {1})".format(
                    self.get_name(), paramstr)


class Gaus95(UncertaintyModel):
    """Gaussian distribution defined by a 95% confidence interval range."""

    def __init__(self, *args):
        if len(args) == 1:
            self.low = -args[0]
            self.high = args[0]
        elif len(args) == 2:
            if args[0] > args[1]:
                msg = "The second argument '{0}' must be bigger than the" \
                      " first '{1}'.".format(args[1], args[0])
                raise SyntaxError(msg)
            self.low = args[0]
            self.high = args[1]
        else:
            msg = "Maximum of 2 arguments supported for the {0}" \
                  " uncertainty model.".format(self.__class__.__name__)
            raise SyntaxError(msg)

    def __repr__(self):
        return "Gaussian uncertainty with 95% confidence interval [{0} .." \
               " {1}]".format(self.low, self.high)

    def get_name(self):
        return "gaus95"

    def get_range(self):
        # Use the 95% confidence interval as the DVH data point range for a
        # Gaussian distribution.
        return (self.high - self.low) * 0.5

    def as_matlab(self):
        paramstr = "{{" + str(self.low) + ", " + str(self.high) + "}}"
        return "struct('uncertainty_model', '{0}', 'params', {1})".format(
                    self.get_name(), paramstr)


class LogNorm(UncertaintyModel):
    """Log-normal distribution defined by parameters mu and sigma."""

    def __init__(self, *args):
        if len(args) == 1:
            self.mu = 0
            self.sigma = args[0]
        elif len(args) == 2:
            self.mu = args[0]
            self.sigma = args[1]
        else:
            msg = "Maximum of 2 arguments supported for the {0}" \
                  " uncertainty model.".format(self.__class__.__name__)
            raise SyntaxError(msg)

    def __repr__(self):
        return "Log normal uncertainty mu = {0}, sigma = {1}".format(self.mu,
                                                                     self.sigma)

    def get_name(self):
        return "lognorm"

    def as_matlab(self):
        paramstr = "{{" + str(self.mu) + ", " + str(self.sigma) + "}}"
        return "struct('uncertainty_model', '{1}', 'params', {2});\n".format(
                    self.get_name(), paramstr)


class LogNorm95(UncertaintyModel):
    """Log-normal distribution defined by a 95% confidence interval range."""

    def __init__(self, low, high):
        if low > high:
            msg = "The second argument '{0}' must be bigger than the" \
                  " first '{1}'.".format(high, low)
            raise SyntaxError(msg)
        if low < 0:
            msg = "The arguments must be greater than zero, but got" \
                  " '{0}'.".format(low)
            raise SyntaxError(msg)
        self.low = low
        self.high = high

    def __repr__(self):
        return "Log normal uncertainty with 95% confidence interval [{0} .." \
               " {1}]".format(self.low, self.high)

    def get_name(self):
        return "lognorm95"

    def as_matlab(self):
        paramstr = "{{" + str(self.low) + ", " + str(self.high) + "}}"
        return "struct('uncertainty_model', '{1}', 'params', {2});\n".format(
                    self.get_name(), paramstr)


class OrganParams(object):
    """Object of model parameters for a particular organ."""

    registery = []

    def __init__(self, name, **kwargs):
        for key in kwargs:
            params = kwargs[key]
            if not isinstance(params, list):
                raise SyntaxError("'{0}' must be assigned a list of uncertanty"
                                  " model objects.".format(key))
            for n, p in enumerate(params):
                if not isinstance(p, UncertaintyModel):
                    msg = "Parameter {0} of '{1}' must be an uncertainty" \
                          " model object.".format(n+1, key)
                    raise SyntaxError(msg)
        self.name = name
        self.models = kwargs
        OrganParams.registery.append(self)


class ConfigParams(object):
    """
    Object to represent configuration parameters loaded from the external
    configuration file.
    """

    def __init__(self, integration_methods, interpolation_methods,
                 dose_binning_uncertainty, volume_ratio_uncertainty,
                 bootstrap_max_samples, bootstrap_sample_mode,
                 organ_name_map, organ_params, organs, models):
        self.integration_methods = integration_methods
        self.interpolation_methods = interpolation_methods
        self.dose_binning_uncertainty = dose_binning_uncertainty
        self.volume_ratio_uncertainty = volume_ratio_uncertainty
        self.bootstrap_max_samples = bootstrap_max_samples
        self.bootstrap_sample_mode = bootstrap_sample_mode
        self.organ_name_map = organ_name_map
        self.organ_parameters = collections.defaultdict(dict)
        for organ in organ_params:
            for model in organ.models:
                params = organ.models[model]
                self.organ_parameters[organ.name][model] = params
        self.organs = organs
        self.models = models

    def __repr__(self):
        # Find the maximum length of organ name strings.
        max_organ_keys_size = max(map(len, self.organ_parameters.keys()))
        # Find the maximum length of model name strings.
        model_keys = []
        for organ in self.organ_parameters:
            model_keys += self.organ_parameters[organ].keys()
        max_model_key_size = max(map(len, model_keys))
        # Prepare the header for the organs/models list.
        organ_output = "{0:>{w1}}  {1:>{w2}}          {2}\n".format(
                        "Organ", "Model", "Parameters",
                        w1 = max_organ_keys_size, w2 = max_model_key_size)
        # Construct the table of organs/models and their parameters.
        for organ in self.organ_parameters:
            for model in self.organ_parameters[organ]:
                params = self.organ_parameters[organ][model]
                if len(params) > 0:
                    pattern = "{0:>{w1}}  {1:>{w2}}  {2:>5}:  {3}\n"
                    for n, param in enumerate(params):
                        organ_output += pattern.format(
                            organ if n == 0 else "", model if n == 0 else "",
                            n+1, param, w1 = max_organ_keys_size,
                            w2 = max_model_key_size)
                else:
                    pattern = "{0:>{w1}}  {1:>{w2}}          (none)\n"
                    organ_output += pattern.format(
                        organ, model, n, w1 = max_organ_keys_size,
                        w2 = max_model_key_size)
        # Construct and return the full string representation:
        return "organs = {0}\n" \
               "models = {1}\n" \
               "integration_methods = {2}\n" \
               "interpolation_methods = {3}\n" \
               "dose_binning_uncertainty = {4}\n" \
               "volume_ratio_uncertainty = {5}\n" \
               "bootstrap_max_samples = {6}\n" \
               "bootstrap_sample_mode = {7}\n" \
               "organ_name_map = {8}\n{9}".format(
                   self.organs, self.models, self.integration_methods,
                   self.interpolation_methods, self.dose_binning_uncertainty,
                   self.volume_ratio_uncertainty, self.bootstrap_max_samples,
                   self.bootstrap_sample_mode, self.organ_name_map,
                   organ_output.rstrip())

    def generate_matlab_params(self):
        """generate_matlab_params() -> string

        Generates a string containing Matlab code to prepare the parameters
        needed for the processPatients.m function.
        """
        script = ""
        # Prepare the organ_name_map variable.
        func = lambda x: "'{0}', '{1}'".format(x, self.organ_name_map[x])
        liststr = "; ".join(map(func, self.organ_name_map))
        script += "organ_name_map = {" + liststr + "};\n"
        # Add the organ list.
        liststr = ", ".join(map(lambda x: "'{0}'".format(x), self.organs))
        script += "organs = {" + liststr + "};\n"
        # Add the model list.
        liststr = ", ".join(map(lambda x: "'{0}'".format(x), self.models))
        script += "models = {" + liststr + "};\n"
        # Prepare the integration and interpolation method code snipets.
        for n, method in enumerate(self.integration_methods):
            tolerance = self.integration_methods[method]
            pattern = "_method{0} = struct('name', '{1}', 'tolerance', {2});\n"
            script += pattern.format(n, method, tolerance)
        liststr = ", ".join(map(lambda x: "_method{0}".format(x),
                            xrange(len(self.integration_methods))))
        integration_methods = "{" + liststr + "}"
        liststr = ", ".join(map(lambda x: "'{0}'".format(x),
                                self.interpolation_methods))
        interpolation_methods = "{" + liststr + "}"
        # Prepare the organ parameter structure.
        if self.organ_parameters:
            script += "_organ_params = struct;\n"
            pnum = 0
            for organ in self.organ_parameters:
                for model in self.organ_parameters[organ]:
                    params = self.organ_parameters[organ][model]
                    paramstr = "{"
                    for n, p in enumerate(params):
                        script += "_model_param{0} = {1};\n".format(
                                                            pnum, p.as_matlab())
                        if n > 0: paramstr += ", "
                        paramstr += "_model_param{0}".format(pnum)
                        pnum += 1
                    paramstr += "}"
                    script += "_organ_params.{0}.{1} = {2};\n".format(
                                                        organ, model, paramstr)
        # Assemble the params structure code.
        if self.dose_binning_uncertainty or self.volume_ratio_uncertainty \
           or self.integration_methods or self.interpolation_methods \
           or self.organ_parameters or self.bootstrap_max_samples \
           or self.bootstrap_sample_mode:
            script += "params = struct;\n"
        else:
            script += "params = {};\n"
        if self.dose_binning_uncertainty:
            script += "params.dose_binning_uncertainty_model = '{0}';\n".format(
                                    self.dose_binning_uncertainty.get_name())
            script += "params.dose_binning_uncertainty = {0};\n".format(
                                    self.dose_binning_uncertainty.get_range())
        if self.volume_ratio_uncertainty:
            script += "params.volume_ratio_uncertainty_model = '{0}';\n".format(
                                    self.volume_ratio_uncertainty.get_name())
            script += "params.volume_ratio_uncertainty = {0};\n".format(
                                    self.volume_ratio_uncertainty.get_range())
        if self.integration_methods:
            script += "params.integration_methods = {0};\n".format(
                                                        integration_methods)
        if self.interpolation_methods:
            script += "params.interpolation_methods = {0};\n".format(
                                                        interpolation_methods)
        if self.bootstrap_max_samples:
            script += "params.bootstrap_max_samples = {0};\n".format(
                                                    self.bootstrap_max_samples)
        if self.bootstrap_sample_mode:
            script += "params.bootstrap_sample_mode = '{0}';\n".format(
                                                    self.bootstrap_sample_mode)
        if self.organ_parameters:
            script += "params.organs = _organ_params;\n"
        return script


def prepare_argument_parser():
    """prepare_argument_parser() -> object

    Prepares the command line argument parser.
    """
    argparser = argparse.ArgumentParser(
        description = """Processes DVH patient data by performing Monte-Carlo
            sampling of the DVH data and response model parameters.
            The results include values such as calculated OEDs.""")
    argparser.add_argument("filelist", metavar = "<file>", nargs = '+',
        help = """DVH input files to process. Must be in Matlab format (.mat)
            that can be converted with the convertDVHtoMatlabFile.sh script from
            text files.""")
    argparser.add_argument("-c", "--config", dest = "configfile",
        default = None, metavar = "<file>", action = SetUnique,
        help = """Provides a configuration file for the organ/model
            parameters that uses python syntax.""")
    argparser.add_argument("-p", "--print", dest = "printconfig",
        default = False, action = "store_true",
        help = """Prints the configuration loaded.""")
    argparser.add_argument("-o", "--organ", dest = "organlist",
        default = [], metavar = "<name>", action = "append",
        help = """A name of an organ to process.""")
    argparser.add_argument("-m", "--model", dest = "modellist",
        default = [], metavar = "<name>", action = "append",
        help = """A response model to calculate.""")
    argparser.add_argument("-O", "--outfile", dest = "outputfile",
        default = "output.mat", metavar = "<file>", action = SetUnique,
        help = """The name of the output file to write.""")
    argparser.add_argument("-N", "--nsamples", dest = "nsamples",
        default = 100, metavar = "<number>", action = SetUnique,
        help = """The number of samples to generate when performing the
            Monte-Carlo sampling.""")
    return argparser


def load_config(filename, organlist, modellist):
    """load_config(string, list, list) -> object

    Load parameter configurations from the given filename. Also need to provide
    default organ and model lists.
    """
    variables = {'Delta': Delta, 'DoubleDelta': DoubleDelta,
                 'Box': Box, 'Box95': Box95, 'Triangle': Triangle,
                 'Triangle95': Triangle95, 'Triangle95mode': Triangle95mode,
                 'Gaus': Gaus, 'Gaus95': Gaus95, 'LogNorm': LogNorm,
                 'LogNorm95': LogNorm95, 'OrganParams': OrganParams}
    try:
        execfile(filename, variables)
        integration_methods = variables['integration_methods']
        if not isinstance(integration_methods, dict):
            raise RunError("'integration_methods' in the configuration file"
                           " '{0}' must be a dictionary.".format(filename))
        interpolation_methods = variables['interpolation_methods']
        if not isinstance(interpolation_methods, list):
            raise RunError("'interpolation_methods' in the configuration file"
                           " '{0}' must be a list.".format(filename))
        dose_binning_uncertainty = variables['dose_binning_uncertainty']
        if not isinstance(dose_binning_uncertainty, UncertaintyModel):
            msg = "'dose_binning_uncertainty' in the configuration file '{0}'" \
                  " must be an uncertainty model object.".format(filename)
            raise RunError(msg)
        volume_ratio_uncertainty = variables['volume_ratio_uncertainty']
        if not isinstance(volume_ratio_uncertainty, UncertaintyModel):
            msg = "'volume_ratio_uncertainty' in the configuration file '{0}'" \
                  " must be an uncertainty model object.".format(filename)
            raise RunError(msg)
        # The following fields are made optional.
        if 'organ_name_map' in variables:
            organ_name_map = variables['organ_name_map']
        else:
            organ_name_map = {}
        if not isinstance(organ_name_map, dict):
            raise RunError("'organ_name_map' in the configuration file"
                           " '{0}' must be a dictionary.".format(filename))
        if 'organs' in variables:
            organs = variables['organs']
        else:
            organs = []
        if not isinstance(organs, list):
            raise RunError("'organs' in the configuration file '{0}' must"
                           " be a list of strings.".format(filename))
        if 'models' in variables:
            models = variables['models']
        else:
            models = []
        if not isinstance(models, list):
            raise RunError("'models' in the configuration file '{0}' must"
                           " be a list of strings.".format(filename))
        if 'bootstrap_max_samples' in variables:
            bootstrap_max_samples = variables['bootstrap_max_samples']
        else:
            bootstrap_max_samples = 6435
        if not isinstance(bootstrap_max_samples, int):
            raise RunError("'bootstrap_max_samples' in the configuration file"
                           " '{0}' must be an integer.".format(filename))
        if 'bootstrap_sample_mode' in variables:
            bootstrap_sample_mode = variables['bootstrap_sample_mode']
        else:
            bootstrap_sample_mode = 'adaptive'
        if not isinstance(bootstrap_sample_mode, str):
            raise RunError("'bootstrap_sample_mode' in the configuration file"
                           " '{0}' must be a string.".format(filename))
    except IOError as e:
        raise RunError(str(e))
    except KeyError as e:
        raise RunError("Missing declaration of variable {0} in the"
                       " configuration file '{1}'.".format(e, filename))
    organs += organlist
    models += modellist
    organs = list(set(organs))  # make unique
    models = list(set(models))  # make unique
    return ConfigParams(integration_methods, interpolation_methods,
                        dose_binning_uncertainty, volume_ratio_uncertainty,
                        bootstrap_max_samples, bootstrap_sample_mode,
                        organ_name_map, OrganParams.registery, organs, models)


def run():
    """run() -> None

    Parse the command line arguments and configuration file and then process
    the patient DVH files.
    """
    argparser = prepare_argument_parser()
    args = argparser.parse_args()
    # Load the configuration file and print it if so requested.
    if args.configfile:
        params = load_config(args.configfile, args.organlist, args.modellist)
    else:
        organs = list(set(args.organlist))
        models = list(set(args.modellist))
        params = ConfigParams({}, [], None, None, None, None, {}, [],
                              organs, models)
    if args.printconfig:
        print("Parameter configuration:")
        if args.configfile:
            print(params)
        else:
            print("(No config given. Using default internal parameters.)")
    # Check that the input files exist.
    for filename in args.filelist:
        if not os.path.exists(filename):
            raise RunError("The file '{0}' does not exist.".format(filename))
    # Prepare the Matlab script to execute, starting with the DVH file list.
    script = ""
    fileliststr = ", ".join(map(lambda x: "'{0}'".format(x), args.filelist))
    script += "dvhfiles = {" + fileliststr + "};\n"
    # Add the parameter script snippet.
    script += params.generate_matlab_params()
    # Add the function call to process the patient input files.
    script += "Results = processPatients(dvhfiles, params, organs, models," \
              " organ_name_map, {0});\n".format(args.nsamples)
    # Add the function calls to merge data if the output file already exists.
    if os.path.exists(args.outputfile):
        script += textwrap.dedent("""\
            try
                Previous = load('{0}', 'Results').Results;
            catch
                error('Missing Results in "{0}".');
            end
            if ~ isstruct(Previous)
                error('Results in "{0}" is not a structure.');
            end
            Results = mergePatientResults(Previous, Results);
            """.format(args.outputfile))
    script += "save('-7', '{0}', 'Results');\n".format(args.outputfile)
    # Invoke octave to execute the Matlab script.
    command_less_script = ['octave', '-q', '--path', '@@MFILE_PATH@@', '--eval']
    command = command_less_script + [script]
    if subprocess.call(command) != 0:
        cmdstr = " ".join(command_less_script + ["'...'"])
        filename = "failed_script.m"
        with open(filename, "w") as scriptfile:
            scriptfile.write(script)
        raise RunError("Execution of command failed: {0}\nWrote the failing"
                       " script that to file '{1}'.".format(cmdstr, filename))


if __name__ == '__main__':
    try:
        run()
    except RunError as e:
        print("ERROR: " + str(e), file = sys.stderr)
