
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <GlobalInc.h>
#include <IM_IO.h>
#include <TempArray.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================

#include <numpy/arrayobject.h>
#define PyArray_SimpleNewFromDataFortran(nd, dims, typenum, data)  	PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                    data, 0, NPY_FARRAY, NULL)
#define PyArray_SimpleNewFromDataC(nd, dims, typenum, data)  	PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                    data, 0, NPY_CARRAY, NULL)

#define generic_array(name,PARAM_TYPE,ARRAY_TYPE, NPY_TYPE)     PyObject* name(to_array<PARAM_TYPE, ARRAY_TYPE> &self){       int ndim = self.naxis();        PyArray_ORDER order = NPY_FORTRANORDER;      npy_intp * dims = new npy_intp[ndim];      for(int i=0; i < ndim; i++) dims[i] = self.axis(i+1);       PyObject *out= PyArray_SimpleNewFromDataFortran(ndim,dims,NPY_TYPE,self.buffer());       delete[] dims;       return out;} 
      
#define generic_image(name,PARAM_TYPE,ARRAY_TYPE, NPY_TYPE)     PyObject* name(to_array<PARAM_TYPE, ARRAY_TYPE> &self){       int ndim = self.naxis();        PyArray_ORDER order = NPY_CORDER;       npy_intp * dims = new npy_intp[ndim];       for(int i=0; i < ndim; i++) dims[i] = self.axis(i+1);       if(ndim >=2 ){ dims[0] = self.nl(); dims[1] = self.nc();}       PyObject *out= PyArray_SimpleNewFromDataC(ndim,dims,NPY_TYPE,self.buffer());       delete[] dims;       return out;} 
      
generic_array(dataf,float,true,NPY_FLOAT32)
generic_array(datad,double,true,NPY_DOUBLE)
generic_image(dataIf,float,false,NPY_FLOAT32)
generic_image(dataId,double,false,NPY_DOUBLE)


namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_alloc_overloads_1_2, alloc, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_alloc_overloads_2_3, alloc, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_alloc_overloads_3_4, alloc, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_alloc_overloads_3_5, alloc, 3, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_alloc_overloads_4_6, alloc, 4, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_reform_overloads_1_3, reform, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_resize_overloads_1_3, resize, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_info_overloads_0_1, info, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_display_overloads_0_1, display, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_sigma_clip_overloads_2_3, sigma_clip, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_true_sigma_clip_overloads_0_1, sigma_clip, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_alloc_overloads_1_2, alloc, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_alloc_overloads_2_3, alloc, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_alloc_overloads_3_4, alloc, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_alloc_overloads_3_5, alloc, 3, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_alloc_overloads_4_6, alloc, 4, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_reform_overloads_1_3, reform, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_resize_overloads_1_3, resize, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_info_overloads_0_1, info, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_display_overloads_0_1, display, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_sigma_clip_overloads_2_3, sigma_clip, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_true_sigma_clip_overloads_0_1, sigma_clip, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_alloc_overloads_1_2, alloc, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_alloc_overloads_2_3, alloc, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_alloc_overloads_3_4, alloc, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_alloc_overloads_3_5, alloc, 3, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_alloc_overloads_4_6, alloc, 4, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_reform_overloads_1_3, reform, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_resize_overloads_1_3, resize, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_info_overloads_0_1, info, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_display_overloads_0_1, display, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_sigma_clip_overloads_2_3, sigma_clip, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_float_false_sigma_clip_overloads_0_1, sigma_clip, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_alloc_overloads_1_2, alloc, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_alloc_overloads_2_3, alloc, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_alloc_overloads_3_4, alloc, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_alloc_overloads_3_5, alloc, 3, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_alloc_overloads_4_6, alloc, 4, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_reform_overloads_1_3, reform, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_resize_overloads_1_3, resize, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_info_overloads_0_1, info, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_display_overloads_0_1, display, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_sigma_clip_overloads_2_3, sigma_clip, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(to_array_double_false_sigma_clip_overloads_0_1, sigma_clip, 0, 1)

BOOST_PYTHON_FUNCTION_OVERLOADS(io_read_ima_float_overloads_2_3, io_read_ima_float, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(io_write_ima_float_overloads_2_4, io_write_ima_float, 2, 4)

}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(isap_great3_tools)
{

    import_array();
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    enum_< Bool >("Bool")
        .value("True", True)
        .value("False", False)
    ;

    class_< to_array<float,true> >("fltarr", init<  >())
        .def(init< int, const char* >())
        .def(init< int, int, const char* >())
        .def(init< int, int, int, const char* >())
        .def(init< int, optional< int, int > >())
        .def(init< const to_array<float,true>& >())
        .def_readwrite("JunkVar", &to_array<float,true>::JunkVar)
        .def("set2ima", &to_array<float,true>::set2ima)
        .def("set2tab", &to_array<float,true>::set2tab)
        .def("free", &to_array<float,true>::free)
        .def("init", (void (to_array<float,true>::*)(const to_array<float,true>&) )&to_array<float,true>::init)
        .def("init", (void (to_array<float,true>::*)() )&to_array<float,true>::init)
        .def("init", (void (to_array<float,true>::*)(float) )&to_array<float,true>::init)
        .def("alloc", (void (to_array<float,true>::*)(int, const char*) )&to_array<float,true>::alloc, to_array_float_true_alloc_overloads_1_2())
        .def("alloc", (void (to_array<float,true>::*)(int, int, const char*) )&to_array<float,true>::alloc, to_array_float_true_alloc_overloads_2_3())
        .def("alloc", (void (to_array<float,true>::*)(int, int, int, const char*) )&to_array<float,true>::alloc, to_array_float_true_alloc_overloads_3_4())
        .def("alloc", (void (to_array<float,true>::*)(float*, int, int, const char*, bool) )&to_array<float,true>::alloc, to_array_float_true_alloc_overloads_3_5())
        .def("alloc", (void (to_array<float,true>::*)(float*, int, int, int, const char*, bool) )&to_array<float,true>::alloc, to_array_float_true_alloc_overloads_4_6())
        .def("reform", &to_array<float,true>::reform, to_array_float_true_reform_overloads_1_3())
        .def("resize", &to_array<float,true>::resize, to_array_float_true_resize_overloads_1_3())
        .def("info", &to_array<float,true>::info, to_array_float_true_info_overloads_0_1())
        .def("display", &to_array<float,true>::display, to_array_float_true_display_overloads_0_1())
        .def("rampgen", &to_array<float,true>::rampgen)
        .def("sup_threshold", &to_array<float,true>::sup_threshold)
        .def("inf_threshold", &to_array<float,true>::inf_threshold)
        .def("n_elem", &to_array<float,true>::n_elem)
        .def("naxis", &to_array<float,true>::naxis)
        .def("axis", &to_array<float,true>::axis)
        .def("nx", &to_array<float,true>::nx)
        .def("ny", &to_array<float,true>::ny)
        .def("nz", &to_array<float,true>::nz)
        .def("get_isbuffer", &to_array<float,true>::get_isbuffer)
        .def("set_isbuffer", &to_array<float,true>::set_isbuffer)
        .def("get_buf", &to_array<float,true>::get_buf)
        .def("get_memalloc", &to_array<float,true>::get_memalloc)
        .def("nc", &to_array<float,true>::nc)
        .def("nl", &to_array<float,true>::nl)
        .def("min", (float (to_array<float,true>::*)() )&to_array<float,true>::min)
        .def("max", (float (to_array<float,true>::*)() )&to_array<float,true>::max)
        .def("maxfabs", (float (to_array<float,true>::*)() )&to_array<float,true>::maxfabs)
        .def("min", (float (to_array<float,true>::*)(int&) )&to_array<float,true>::min)
        .def("max", (float (to_array<float,true>::*)(int&) )&to_array<float,true>::max)
        .def("maxfabs", (float (to_array<float,true>::*)(int&) )&to_array<float,true>::maxfabs)
        .def("total", &to_array<float,true>::total)
        .def("energy", &to_array<float,true>::energy)
        .def("sigma", &to_array<float,true>::sigma)
        .def("mean", &to_array<float,true>::mean)
        .def("sigma_clip", (void (to_array<float,true>::*)(float&, float&, int) const)&to_array<float,true>::sigma_clip, to_array_float_true_sigma_clip_overloads_2_3())
        .def("sigma_clip", (float (to_array<float,true>::*)(int) const)&to_array<float,true>::sigma_clip, to_array_float_true_sigma_clip_overloads_0_1())
        .def( self += self )
        .def( self *= self )
        .def( self *= other< double >() )
        .def( self -= self )
        .def( self /= self )
        .def( self ^ other< double >() )
        .add_property("data",&dataf)
    ;

    class_< to_array<double,true> >("dblarr", init<  >())
        .def(init< int, const char* >())
        .def(init< int, int, const char* >())
        .def(init< int, int, int, const char* >())
        .def(init< int, optional< int, int > >())
        .def(init< const to_array<double,true>& >())
        .def_readwrite("JunkVar", &to_array<double,true>::JunkVar)
        .def("set2ima", &to_array<double,true>::set2ima)
        .def("set2tab", &to_array<double,true>::set2tab)
        .def("free", &to_array<double,true>::free)
        .def("init", (void (to_array<double,true>::*)(const to_array<double,true>&) )&to_array<double,true>::init)
        .def("init", (void (to_array<double,true>::*)() )&to_array<double,true>::init)
        .def("init", (void (to_array<double,true>::*)(double) )&to_array<double,true>::init)
        .def("alloc", (void (to_array<double,true>::*)(int, const char*) )&to_array<double,true>::alloc, to_array_double_true_alloc_overloads_1_2())
        .def("alloc", (void (to_array<double,true>::*)(int, int, const char*) )&to_array<double,true>::alloc, to_array_double_true_alloc_overloads_2_3())
        .def("alloc", (void (to_array<double,true>::*)(int, int, int, const char*) )&to_array<double,true>::alloc, to_array_double_true_alloc_overloads_3_4())
        .def("alloc", (void (to_array<double,true>::*)(double*, int, int, const char*, bool) )&to_array<double,true>::alloc, to_array_double_true_alloc_overloads_3_5())
        .def("alloc", (void (to_array<double,true>::*)(double*, int, int, int, const char*, bool) )&to_array<double,true>::alloc, to_array_double_true_alloc_overloads_4_6())
        .def("reform", &to_array<double,true>::reform, to_array_double_true_reform_overloads_1_3())
        .def("resize", &to_array<double,true>::resize, to_array_double_true_resize_overloads_1_3())
        .def("info", &to_array<double,true>::info, to_array_double_true_info_overloads_0_1())
        .def("display", &to_array<double,true>::display, to_array_double_true_display_overloads_0_1())
        .def("rampgen", &to_array<double,true>::rampgen)
        .def("sup_threshold", &to_array<double,true>::sup_threshold)
        .def("inf_threshold", &to_array<double,true>::inf_threshold)
        .def("n_elem", &to_array<double,true>::n_elem)
        .def("naxis", &to_array<double,true>::naxis)
        .def("axis", &to_array<double,true>::axis)
        .def("nx", &to_array<double,true>::nx)
        .def("ny", &to_array<double,true>::ny)
        .def("nz", &to_array<double,true>::nz)
        .def("get_isbuffer", &to_array<double,true>::get_isbuffer)
        .def("set_isbuffer", &to_array<double,true>::set_isbuffer)
        .def("get_buf", &to_array<double,true>::get_buf)
        .def("get_memalloc", &to_array<double,true>::get_memalloc)
        .def("nc", &to_array<double,true>::nc)
        .def("nl", &to_array<double,true>::nl)
        .def("min", (double (to_array<double,true>::*)() )&to_array<double,true>::min)
        .def("max", (double (to_array<double,true>::*)() )&to_array<double,true>::max)
        .def("maxfabs", (double (to_array<double,true>::*)() )&to_array<double,true>::maxfabs)
        .def("min", (double (to_array<double,true>::*)(int&) )&to_array<double,true>::min)
        .def("max", (double (to_array<double,true>::*)(int&) )&to_array<double,true>::max)
        .def("maxfabs", (double (to_array<double,true>::*)(int&) )&to_array<double,true>::maxfabs)
        .def("total", &to_array<double,true>::total)
        .def("energy", &to_array<double,true>::energy)
        .def("sigma", &to_array<double,true>::sigma)
        .def("mean", &to_array<double,true>::mean)
        .def("sigma_clip", (void (to_array<double,true>::*)(float&, float&, int) const)&to_array<double,true>::sigma_clip, to_array_double_true_sigma_clip_overloads_2_3())
        .def("sigma_clip", (float (to_array<double,true>::*)(int) const)&to_array<double,true>::sigma_clip, to_array_double_true_sigma_clip_overloads_0_1())
        .def( self += self )
        .def( self *= self )
        .def( self *= other< double >() )
        .def( self -= self )
        .def( self /= self )
        .def( self ^ other< double >() )
        .add_property("data",&datad)
    ;

    class_< to_array<float,false> >("Iflt", init<  >())
        .def(init< int, const char* >())
        .def(init< int, int, const char* >())
        .def(init< int, int, int, const char* >())
        .def(init< int, optional< int, int > >())
        .def(init< const to_array<float,false>& >())
        .def_readwrite("JunkVar", &to_array<float,false>::JunkVar)
        .def("set2ima", &to_array<float,false>::set2ima)
        .def("set2tab", &to_array<float,false>::set2tab)
        .def("free", &to_array<float,false>::free)
        .def("init", (void (to_array<float,false>::*)(const to_array<float,false>&) )&to_array<float,false>::init)
        .def("init", (void (to_array<float,false>::*)() )&to_array<float,false>::init)
        .def("init", (void (to_array<float,false>::*)(float) )&to_array<float,false>::init)
        .def("alloc", (void (to_array<float,false>::*)(int, const char*) )&to_array<float,false>::alloc, to_array_float_false_alloc_overloads_1_2())
        .def("alloc", (void (to_array<float,false>::*)(int, int, const char*) )&to_array<float,false>::alloc, to_array_float_false_alloc_overloads_2_3())
        .def("alloc", (void (to_array<float,false>::*)(int, int, int, const char*) )&to_array<float,false>::alloc, to_array_float_false_alloc_overloads_3_4())
        .def("alloc", (void (to_array<float,false>::*)(float*, int, int, const char*, bool) )&to_array<float,false>::alloc, to_array_float_false_alloc_overloads_3_5())
        .def("alloc", (void (to_array<float,false>::*)(float*, int, int, int, const char*, bool) )&to_array<float,false>::alloc, to_array_float_false_alloc_overloads_4_6())
        .def("reform", &to_array<float,false>::reform, to_array_float_false_reform_overloads_1_3())
        .def("resize", &to_array<float,false>::resize, to_array_float_false_resize_overloads_1_3())
        .def("info", &to_array<float,false>::info, to_array_float_false_info_overloads_0_1())
        .def("display", &to_array<float,false>::display, to_array_float_false_display_overloads_0_1())
        .def("rampgen", &to_array<float,false>::rampgen)
        .def("sup_threshold", &to_array<float,false>::sup_threshold)
        .def("inf_threshold", &to_array<float,false>::inf_threshold)
        .def("n_elem", &to_array<float,false>::n_elem)
        .def("naxis", &to_array<float,false>::naxis)
        .def("axis", &to_array<float,false>::axis)
        .def("nx", &to_array<float,false>::nx)
        .def("ny", &to_array<float,false>::ny)
        .def("nz", &to_array<float,false>::nz)
        .def("get_isbuffer", &to_array<float,false>::get_isbuffer)
        .def("set_isbuffer", &to_array<float,false>::set_isbuffer)
        .def("get_buf", &to_array<float,false>::get_buf)
        .def("get_memalloc", &to_array<float,false>::get_memalloc)
        .def("nc", &to_array<float,false>::nc)
        .def("nl", &to_array<float,false>::nl)
        .def("min", (float (to_array<float,false>::*)() )&to_array<float,false>::min)
        .def("max", (float (to_array<float,false>::*)() )&to_array<float,false>::max)
        .def("maxfabs", (float (to_array<float,false>::*)() )&to_array<float,false>::maxfabs)
        .def("min", (float (to_array<float,false>::*)(int&) )&to_array<float,false>::min)
        .def("max", (float (to_array<float,false>::*)(int&) )&to_array<float,false>::max)
        .def("maxfabs", (float (to_array<float,false>::*)(int&) )&to_array<float,false>::maxfabs)
        .def("total", &to_array<float,false>::total)
        .def("energy", &to_array<float,false>::energy)
        .def("sigma", &to_array<float,false>::sigma)
        .def("mean", &to_array<float,false>::mean)
        .def("sigma_clip", (void (to_array<float,false>::*)(float&, float&, int) const)&to_array<float,false>::sigma_clip, to_array_float_false_sigma_clip_overloads_2_3())
        .def("sigma_clip", (float (to_array<float,false>::*)(int) const)&to_array<float,false>::sigma_clip, to_array_float_false_sigma_clip_overloads_0_1())
        .def( self += self )
        .def( self *= self )
        .def( self *= other< double >() )
        .def( self -= self )
        .def( self /= self )
        .def( self ^ other< double >() )
        .add_property("data",&dataIf)
    ;

    class_< to_array<double,false> >("Idbl", init<  >())
        .def(init< int, const char* >())
        .def(init< int, int, const char* >())
        .def(init< int, int, int, const char* >())
        .def(init< int, optional< int, int > >())
        .def(init< const to_array<double,false>& >())
        .def_readwrite("JunkVar", &to_array<double,false>::JunkVar)
        .def("set2ima", &to_array<double,false>::set2ima)
        .def("set2tab", &to_array<double,false>::set2tab)
        .def("free", &to_array<double,false>::free)
        .def("init", (void (to_array<double,false>::*)(const to_array<double,false>&) )&to_array<double,false>::init)
        .def("init", (void (to_array<double,false>::*)() )&to_array<double,false>::init)
        .def("init", (void (to_array<double,false>::*)(double) )&to_array<double,false>::init)
        .def("alloc", (void (to_array<double,false>::*)(int, const char*) )&to_array<double,false>::alloc, to_array_double_false_alloc_overloads_1_2())
        .def("alloc", (void (to_array<double,false>::*)(int, int, const char*) )&to_array<double,false>::alloc, to_array_double_false_alloc_overloads_2_3())
        .def("alloc", (void (to_array<double,false>::*)(int, int, int, const char*) )&to_array<double,false>::alloc, to_array_double_false_alloc_overloads_3_4())
        .def("alloc", (void (to_array<double,false>::*)(double*, int, int, const char*, bool) )&to_array<double,false>::alloc, to_array_double_false_alloc_overloads_3_5())
        .def("alloc", (void (to_array<double,false>::*)(double*, int, int, int, const char*, bool) )&to_array<double,false>::alloc, to_array_double_false_alloc_overloads_4_6())
        .def("reform", &to_array<double,false>::reform, to_array_double_false_reform_overloads_1_3())
        .def("resize", &to_array<double,false>::resize, to_array_double_false_resize_overloads_1_3())
        .def("info", &to_array<double,false>::info, to_array_double_false_info_overloads_0_1())
        .def("display", &to_array<double,false>::display, to_array_double_false_display_overloads_0_1())
        .def("rampgen", &to_array<double,false>::rampgen)
        .def("sup_threshold", &to_array<double,false>::sup_threshold)
        .def("inf_threshold", &to_array<double,false>::inf_threshold)
        .def("n_elem", &to_array<double,false>::n_elem)
        .def("naxis", &to_array<double,false>::naxis)
        .def("axis", &to_array<double,false>::axis)
        .def("nx", &to_array<double,false>::nx)
        .def("ny", &to_array<double,false>::ny)
        .def("nz", &to_array<double,false>::nz)
        .def("get_isbuffer", &to_array<double,false>::get_isbuffer)
        .def("set_isbuffer", &to_array<double,false>::set_isbuffer)
        .def("get_buf", &to_array<double,false>::get_buf)
        .def("get_memalloc", &to_array<double,false>::get_memalloc)
        .def("nc", &to_array<double,false>::nc)
        .def("nl", &to_array<double,false>::nl)
        .def("min", (double (to_array<double,false>::*)() )&to_array<double,false>::min)
        .def("max", (double (to_array<double,false>::*)() )&to_array<double,false>::max)
        .def("maxfabs", (double (to_array<double,false>::*)() )&to_array<double,false>::maxfabs)
        .def("min", (double (to_array<double,false>::*)(int&) )&to_array<double,false>::min)
        .def("max", (double (to_array<double,false>::*)(int&) )&to_array<double,false>::max)
        .def("maxfabs", (double (to_array<double,false>::*)(int&) )&to_array<double,false>::maxfabs)
        .def("total", &to_array<double,false>::total)
        .def("energy", &to_array<double,false>::energy)
        .def("sigma", &to_array<double,false>::sigma)
        .def("mean", &to_array<double,false>::mean)
        .def("sigma_clip", (void (to_array<double,false>::*)(float&, float&, int) const)&to_array<double,false>::sigma_clip, to_array_double_false_sigma_clip_overloads_2_3())
        .def("sigma_clip", (float (to_array<double,false>::*)(int) const)&to_array<double,false>::sigma_clip, to_array_double_false_sigma_clip_overloads_0_1())
        .def( self += self )
        .def( self *= self )
        .def( self *= other< double >() )
        .def( self -= self )
        .def( self /= self )
        .def( self ^ other< double >() )
        .add_property("data",&dataId)
    ;

    def("fits_read_fltarr", (void (*)(char*, to_array<float,true>&, fitsstruct*, int))&fits_read_fltarr);
    def("fits_read_fltarr", (void (*)(char*, to_array<float,true>&, fitsstruct*))&fits_read_fltarr);
    def("fits_read_fltarr", (void (*)(char*, to_array<float,true>&))&fits_read_fltarr);
    def("fits_read_dblarr", &fits_read_dblarr);
    def("fits_write_fltarr", (void (*)(char*, to_array<float,true>&))&fits_write_fltarr);
    def("fits_write_fltarr", (void (*)(char*, to_array<float,true>&, fitsstruct*))&fits_write_fltarr);
    def("fits_write_dblarr", &fits_write_dblarr);
    def("io_read_ima_float", &io_read_ima_float, io_read_ima_float_overloads_2_3());
    def("io_write_ima_float", &io_write_ima_float, io_write_ima_float_overloads_2_4());
}

