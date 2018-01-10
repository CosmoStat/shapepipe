"""
Created on Wed May 17 10:22:28 2017

@package pse.pse_header header class
@file pse_header.py
Header class
"""

# -- Python imports
import re


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
    

    def param_value(self, param=None):
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


    def _get_value(self, param=None):
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


    def _operate(self, param=None, param_split=None):
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

        try:
            if len(re.split(op,param))==1:
                return self._get_value(param)
            param_add = re.split('\+',param)
            if len(param_add) > 2:
                tmp=0
                for i in param_add:
                    tmp=tmp+self._operate(i,param_split)
                return tmp
            elif len(param_add) == 2:
                if (param_add[0] in param_split) & (param_add[1] in param_split):
                    return self._get_value(param_add[0])+self._get_value(param_add[1])
                else:
                    if param_add[0] in param_split:
                        return self._get_value(param_add[0])+self._operate(param_add[1],param_split)
                    elif param_add[1] in param_split:
                        return self._operate(param_add[0],param_split)+self._get_value(param_add[1])
                    else:
                        return self._operate(param_add[0],param_split)+self._operate(param_add[1],param_split)
            else:
                param_sub = re.split('\-',param)
                if len(param_sub) > 2:
                    tmp='init'
                    for i in param_sub:
                        if tmp=='init':
                            tmp=self._operate(i,param_split)
                        else:
                            tmp=tmp-self._operate(i,param_split)
                    return tmp
                elif len(param_sub) == 2:
                    if (param_sub[0] in param_split) & (param_sub[1] in param_split):
                        return self._get_value(param_sub[0])-self._get_value(param_sub[1])
                    else:
                        if param_sub[0] in param_split:
                            return self._get_value(param_sub[0])-self._operate(param_sub[1],param_split)
                        elif param_sub[1] in param_split:
                            return self._operate(param_sub[0],param_split)-self._get_value(param_sub[1])
                        else:
                            return self._operate(param_sub[0],param_split)-self._operate(param_sub[1],param_split)
                else:
                    param_mul = re.split('\*',param)
                    if len(param_mul) > 2:
                        tmp=1
                        for i in param_mul:
                            tmp=tmp*self._operate(i,param_split)
                        return tmp
                    elif len(param_mul) == 2:
                        if (param_mul[0] in param_split) & (param_mul[1] in param_split):
                            return self._get_value(param_mul[0])*self._get_value(param_mul[1])
                        else:
                            if param_mul[0] in param_split:
                                return self._get_value(param_mul[0])*self._operate(param_mul[1],param_split)
                            elif param_mul[1] in param_split:
                                return self._operate(param_mul[0],param_split)*self._get_value(param_mul[1])
                            else:
                                return self._operate(param_mul[0],param_split)*self._operate(param_mul[1],param_split)
                    else:
                        param_div = re.split('\/',param)
                        if len(param_div) > 2:
                            tmp='init'
                            for i in param_div:
                                if tmp == 'init':
                                    tmp=self._operate(i,param_split)
                                else:
                                    tmp=tmp/self._operate(i,param_split)
                            return tmp
                        elif len(param_div) == 2:
                            if (param_div[0] in param_split) & (param_div[1] in param_split):
                                return self._get_value(param_div[0])/self._get_value(param_div[1])
                            else:
                                if param_div[0] in param_split:
                                    return self._get_value(param_div[0])/self._operate(param_div[1],param_split)
                                elif param_div[1] in param_split:
                                    return self._operate(param_div[0],param_split)/self._get_value(param_div[1])
                                else:
                                    return self._operate(param_div[0],param_split)/self._operate(param_div[1],param_split)
        except:
            Exception("error occurred")


# -- EOF header_handling.py