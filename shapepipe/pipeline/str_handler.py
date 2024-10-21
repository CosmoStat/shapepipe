"""STRING HANDLER.

This module defines a class for handling string operations.

:Author: Axel Guinot

"""

import itertools
import operator
import re

import numpy as np
from modopt.math.stats import sigma_mad


class StrInterpreter(object):
    """String Interpreter Class.

    Class to handle operation/comparison in a string.

    Parameters
    ----------
    string : str
        String to interpret
    catralog : dict, recarray or str
        If type(catalogue) == str : assume a SExtractor fits catalogue and
        read parameter in it else : assume the catalogue is already open and
        look into it
    make_compare : bool
        If true assume a comparison in the string
    mask_dict : dict
        Dictionary containing mask usable for the operation
    autorun : bool
        If true return directly the result

    """

    def __init__(self, string, catalogue, make_compare=False, mask_dict=None):

        if type(string) is not str:
            raise ValueError("string has to be str type")
        else:
            self._string = string

        if catalogue is not None:
            if type(catalogue) is str:
                file = FITSCatalogue(catalogue, SEx_catalogue=True)
                file.open()
                self._cat = file.get_data()
                file.close()
            else:
                self._cat = catalogue
        else:
            raise ValueError("catalogue not provided")

        self._make_compare = make_compare

        if mask_dict is not None:
            self._mask = mask_dict

        self._init_stat_function()
        self._comp_dict = {
            "<": operator.lt,
            ">": operator.gt,
            "<=": operator.le,
            ">=": operator.ge,
            "==": operator.eq,
            "!=": operator.ne,
        }

        self.result = self.interpret(self._string, self._make_compare)

    def interpret(
        self,
        string,
        make_compare=False,
        make_func=True,
        make_operate=True,
    ):
        """Interpret.

        This function handles the different possible operations

        Parameters
        ----------
        str:  string
            string to interpret
        make_compare : bool
            If true look for a comparison
        make_func : bool
            If true look for a function
        make_operate : bool
            If true look for an operation

        Returns
        -------
        array or float
            Result of the current operation.

        Notes
        -----
        This is a recursive function.

        """
        if make_compare:
            result = self._compare(string)
        else:
            if make_operate:
                string_split = re.split(r"\*|\/|\-|\+\s*(?![^()]*\))", string)
                result = self._operate(string, string_split)
            else:
                if make_func:
                    result = self._apply_func(string)
                else:
                    result = self._get_value(string)

        return result

    def _compare(self, string):
        """Handle Comparison in a String.

        This function transform condition in a string to real condition.

        Parameters
        ----------
        string : str
            strind containing the comparison.

        """
        comp = "<|>|<=|>=|==|!="

        if len(re.split(comp, string)) != 2:
            raise Exception(
                "Only one comparison in [<, >, <=, >=, ==, !=] per line"
            )

        for operator in ["<=", ">=", "<", ">", "==", "!="]:
            terms = re.split(operator, string)
            if len(terms) == 2:
                self._make_compare = False
                first = self.interpret(terms[0], self._make_compare)
                second = self.interpret(terms[1], self._make_compare)

                return self._comp_dict[operator](first, second)

    def _apply_func(self, string):
        """Parse Input String for Function Name and Apply Function.

        Parameters
        ----------
        str: string
            input string

        Returns
        -------
        float
            result of the function

        """
        str_split = re.split(r"\(|\)", string)

        if len(str_split) == 1:
            return self.interpret(
                str_split[0],
                self._make_compare,
                make_func=False,
                make_operate=False,
            )
        elif len(str_split) == 3:
            str_split_2 = re.split(",", str_split[1])
            if len(str_split_2) > 1:
                param = [
                    self.interpret(
                        char,
                        self._make_compare,
                        make_func=False,
                        make_operate=True,
                    )
                    for char in str_split_2
                ]

                # Evaluate statistical function, raise error if failure
                # occurs during computation
                try:
                    res = self._stat_func[str_split[0]](*param)
                except Exception:
                    raise
                return res

            else:
                if str_split[0] not in self._stat_func:
                    raise KeyError(
                        f"Invalid function '{str_split[0]}' in expression "
                        + f"'{string}'"
                    )
                return self._stat_func[str_split[0]](
                    self.interpret(
                        str_split[1],
                        self._make_compare,
                        make_func=False,
                        make_operate=True,
                    )
                )
        else:
            raise Exception(
                "Only one function can be applied. Problem with the "
                + f"term: {string}"
            )

    def _init_stat_function(self):
        """Initialise Available Stat Functions.

        Create a dictionary containing the functions.

        """
        self._stat_func = {}
        self._stat_func["mean"] = np.mean
        self._stat_func["median"] = np.median
        self._stat_func["mode"] = self._mode
        self._stat_func["sqrt"] = np.sqrt
        self._stat_func["pow"] = pow
        self._stat_func["log"] = np.log
        self._stat_func["log10"] = np.log10
        self._stat_func["exp"] = np.exp
        self._stat_func["std"] = np.std
        self._stat_func["var"] = np.var
        self._stat_func["sigma_mad"] = self._sigma_mad
        self._stat_func["len"] = len
        self._stat_func["min"] = min
        self._stat_func["max"] = max
        self._stat_func["homogen"] = self._test_homogeneity

    def _mean(self, input):
        """Get Mean.

        Compute the mean of a distribution.

        Parameters
        ----------
        input : numpy.ndarray
            Numpy array containing the data.

        Returns
        -------
        float
            mean, if input array has size > 0; -1, otherwise

        """
        cat_size = len(input)
        if cat_size == 0:
            return -1
        else:
            return np.mean()

    def _mode(self, input, eps=0.001, iter_max=1000):
        """Get Mode.

        Compute the mode, the most frequent value of a continuous distribution.

        Parameters
        ----------
        input : numpy.ndarray
            Numpy array containing the data.
        eps : float, optional
            Accuracy to achieve (default is 0.001)
        iter_max : int, optional
            Maximum number of iterations

        Returns
        -------
        float
            mode, if input array has 10 or more elements;
            median, if input array has >0 and <10 elements;
            -1, if input array has 0 elements

        """
        cat_size = len(input)
        if cat_size > 100:
            bins = int(float(cat_size) / 10.0)
        elif cat_size >= 20:
            bins = int(float(cat_size) / 5.0)
        elif cat_size > 0:
            return np.median(input)
        else:
            return -1

        data = input
        diff = eps + 1.0

        iteration = 0
        while diff > eps:
            hist = np.histogram(data, bins)
            if hist[0].max() == 1:
                break

            b_min = hist[1][hist[0].argmax()]
            b_max = hist[1][hist[0].argmax() + 1]

            diff = b_max - b_min

            data = data[(data > b_min) & (data < b_max)]

            if iteration == iter_max:
                break
            iteration += 1

        if iteration == iter_max:
            raise ValueError("Mode computation failed")
        else:
            mode = (b_min + b_max) / 2.0
            return mode

    def _sigma_mad(self, input):
        """Get Mean Absolute Deviation.

        Compute median absolute deviation (MAD).

        Parameters
        ----------
        input : numpy.ndarray
            input data

        Returns
        -------
        float
            MAD, if input size > 0;
            -1 if input size is 0

        """
        if len(input) == 0:
            return -1
        else:
            return sigma_mad(input)

    def _test_homogeneity(self, *args):
        """Test Homogeneity.

        Test homogeneity on 1D or 2D space.

        Parameters
        ----------
        param1 : numpy.ndarray
            Array on which the homogeneity test is performed
        param2 : numpy.ndarray, optional
            Array on which the homogeneity test is performed
        nb_cells : int
            Number of cells in the space. (note : has to be a square number)

        Returns
        -------
        float
            Percentage of inhomogeneity compared to worse possible case
            (based on the standard deviation)

        """
        if len(args) == 2:
            n_param = 1
            param = [args[0]]
            n_cells = args[1]
        elif len(args) == 3:
            n_param = 2
            param = [args[0], args[1]]
            n_cells = args[2]
        else:
            raise ValueError(
                "Inputs should be param_1, param_2 [optional], n_cells"
            )

        if n_param == 2:
            if len(param[0]) != len(param[1]):
                raise ValueError(
                    "Both param_1 and param_2 must have the same "
                    + f"length : {len(param[0])}, {len(param[1])}"
                )

        if np.sqrt(n_cells) % 1 != 0:
            raise ValueError("N_cells must be a square number")

        n_tot = len(param[0])
        homo_ratio = float(n_tot) / float(n_cells)

        param_min = []
        param_max = []
        for idx in param:
            step = (np.max(idx) - np.min(idx)) / pow(
                n_cells, 1.0 / float(n_param)
            )
            param_min.append(
                [val for val in np.arange(np.min(idx), np.max(idx), step)]
            )
            param_max.append(
                [
                    val
                    for val in np.arange(
                        np.min(idx) + step, np.max(idx) + step, step
                    )
                ]
            )

        if n_param == 1:
            n_obj = np.asarray(
                [
                    float(
                        len(
                            np.where(
                                (param[0] >= param_min[0][idx])
                                & (param[0] <= param_max[0][idx])
                            )[0]
                        )
                    )
                    for idx in range(int(n_cells))
                ]
            )
        elif n_param == 2:
            it = itertools.product(range(int(np.sqrt(n_cells))), repeat=2)
            n_obj = np.asarray(
                [
                    float(
                        len(
                            np.where(
                                (param[0] >= param_min[0][idx_i])
                                & (param[0] <= param_max[0][idx_i])
                                & (param[1] >= param_min[1][idx_j])
                                & (param[1] <= param_max[1][idx_j])
                            )[0]
                        )
                    )
                    for idx_i, idx_j in it
                ]
            )

        actual_std = np.std(n_obj / homo_ratio)

        worse_std = np.zeros((int(n_cells), 1))
        worse_std[0] = n_tot / homo_ratio
        worse_std = np.std(worse_std)

        return actual_std / worse_std * 100.0

    def _operate(self, string, string_split):
        """Handle Operation in a String.

        Make operation between catalogue's parameters and/or numbers.

        Parameters
        ----------
        string : str
            Parameter or linear combination of parameters.
        string_split : str
            String split option

        Returns
        -------
        float
            Result of the operation

        Notes
        -----
        It's used as a recursive function

        """
        op = r"\*|\/|\-|\+\s*(?![^()]*\))"
        if string is None:
            raise ValueError("Parameter not specified")
        if string_split is None:
            raise ValueError("Parameters splited not specified")

        if len(re.split(op, string)) == 1:
            return self.interpret(string, make_operate=False)

        tmp = self._string_op_func(
            re.split(r"\+\s*(?![^()]*\))", string),
            string_split,
            operator.add,
            0,
        )
        if not np.isscalar(tmp) or tmp != "pass":
            return tmp
        else:
            tmp = self._string_op_func(
                re.split(r"\-\s*(?![^()]*\))", string),
                string_split,
                operator.sub,
                "init",
            )
            if not np.isscalar(tmp) or tmp != "pass":
                return tmp
            else:
                tmp = self._string_op_func(
                    re.split(r"\*\s*(?![^()]*\))", string),
                    string_split,
                    operator.mul,
                    1,
                )
                if not np.isscalar(tmp) or tmp != "pass":
                    return tmp
                else:
                    return self._string_op_func(
                        re.split(r"\/\s*(?![^()]*\))", string),
                        string_split,
                        operator.truediv,
                        "init",
                    )

    def _string_op_func(self, string_op, string_split, op, tmp):
        r"""Perform a Specified Operation.

        This function handle the posible operation between parameters.

        Parameters
        ----------
        string_op : list
            List of parameters to operate.
        string_split : list
            The different parameter splitted using
            '\*|\/|\-|\+\s*(?![^()]*\))' as delimiter.
        op : func
            The kind of operation provided as an operator function
            (Example : operator.sub).
        tmp : str or float
            Temporary result of the global operation or value to
            initiate operation.

        Returns
        -------
        float or str
            Result of the operation or 'pass' if there are remaining operations

        """
        if len(string_op) > 2:
            for operator in string_op:
                if tmp == "init":
                    tmp = self._operate(operator, string_split)
                else:
                    tmp = op(tmp, self._operate(operator, string_split))
            return tmp
        elif len(string_op) == 2:
            if string_op[0] in string_split:
                first = self.interpret(string_op[0], make_operate=False)
            else:
                first = self._operate(string_op[0], string_split)
            if string_op[1] in string_split:
                second = self.interpret(string_op[1], make_operate=False)
            else:
                second = self._operate(string_op[1], string_split)
            return op(first, second)
        else:
            return "pass"

    def _get_value(self, string):
        """Get Value.

        Return the value of the corresponding parameter. Or return a float
        with a number as parameter.

        Parameters
        ----------
        string : str
            Parameter of the catalogue.

        Returns
        -------
        float
            Value of the parameter. Or float

        Notes
        -----
        You can't perform operations here!

        """
        if string is None:
            raise ValueError("Parameter not specified")

        try:
            string_value = float(string)
            return string_value
        except Exception:
            str_split = re.split(r"\{|\}", string)
            if len(str_split) == 1:
                try:
                    return self._cat[string]
                except Exception:
                    raise ValueError(
                        "String has to be a float or a catalogue parameter. "
                        + f"{string} not found"
                    )
            if len(str_split) == 3:
                if str_split[1] in self._mask.keys():
                    try:
                        return self._cat[str_split[0]][self._mask[str_split[1]]]
                    except Exception:
                        raise ValueError(
                            "String has to be a catalogue parameter. "
                            + f"{str_split[0]} not found"
                        )
                else:
                    raise ValueError(
                        f"Mask has to be provided. {str_split[1]} not "
                        + "found in mask"
                    )
