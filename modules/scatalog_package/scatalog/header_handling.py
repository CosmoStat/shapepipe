# -- Python imports
import re
import operator


#---------------------------------------------------------------------------------------------------
class Header(object):
    """!
        Tools to make operation with parameters in the header
    """

    #-----------------------------------
    def __init__(self, header):

        self._header = header


    #######################
    #   Public methods  #
    #######################


    def param_value(self, param):
        """!
            Return the value of a parameter from the header
            @param param parameter of the header in str type (see NOTE)
            @return value of the parameter (see NOTE)

            NOTE : you can give a number as parameter it will return this number as float.
                   Also, you can enter a linear combinaition of parameters with numbers in it.
            EXAMPLE : 'GAIN+NEXP*TEXP+10'

            NOTE 2: if the parameter is 'None' link to a DEFAULT value !
        """

        if self._header is None:
            raise ValueError("Header not specified")
        if param is None:
            raise ValueError("Parameter not specified")

        param=re.sub(' ','',param)
        param_split=re.split('\*|\/|\-|\+',param)
        if len(param_split)==1:
            return self._get_value(param)
        else:
            return self._operate(param,param_split)


    #########################
    #   Private methods   #
    #########################


    def _get_value(self, param):
        """!
            Return the value of the corresponding parameter. Or return a float with a number as parameter.
            @param param parameter of the header in str type
            @return value of the parameter. Or float

            NOTE : you can't make operation here !
        """

        if param is None:
            raise ValueError("Parameter not specified")

        try:
            param_value=float(param)
            return param_value
        except:
            try:
                return self._header[param]
            except:
                raise ValueError('param has to be a float or an header parameter')


    def _operate(self, param, param_split):
        """!
            Make operation between header's parameters and/or numbers
            @param param parameter or linear combination of parameters
            @param param_split the different parameter splitted using '\*|\/|\-|\+' as delimiter
            @return result of the operation

            NOTE : it's used as a recursive function
        """

        op='\*|\/|\-|\+'
        if param is None:
            raise ValueError("Parameter not specified")
        if param_split is None:
            raise ValueError("Parameters splited not specified")

        if len(re.split(op,param))==1:
            return self._get_value(param)

        tmp = self._param_op_func(re.split('\+',param), param_split, operator.add, 0)
        if tmp != 'pass':
            return tmp
        else:
            tmp = self._param_op_func(re.split('\-',param), param_split, operator.sub, 'init')
            if tmp != 'pass':
                return tmp
            else:
                tmp = self._param_op_func(re.split('\*',param), param_split, operator.mul, 1)
                if tmp != 'pass':
                    return tmp
                else:
                    return self._param_op_func(re.split('\/',param), param_split, operator.div, 'init')


    def _param_op_func(self, param_op, param_split, op, tmp):
        """!
            This function handle the posible operation between parameters
            @param param_op list of parameters to operate
            @param param_split the different parameter splitted using '\*|\/|\-|\+' as delimiter
            @param op the kind of operation provide as an operator function (Example : operator.sub)
            @param tmp temporary result of the global operation
            @return result of the operation
        """

        if len(param_op) > 2:
            for i in param_op:
                if tmp == 'init':
                    tmp = self._operate(i, param_split)
                else:
                    tmp = op(tmp, self._operate(i, param_split))
            return tmp
        elif len(param_op) == 2:
            if param_op[0] in param_split:
                first = self._get_value(param_op[0])
            else:
                first = self._operate(param_op[0], param_split)
            if param_op[1] in param_split:
                second = self._get_value(param_op[1])
            else:
                second = self._operate(param_op[1], param_split)
            return op(first, second)
        else:
            return 'pass'

# -- EOF header_handling.py
