
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <Border.h>
#include <Filter.h>
#include <IM_Noise.h>
#include <IM_Obj.h>
#include <MR_Filter.h>
#include <MR_NoiseModel.h>
#include <MR_Obj.h>
#include <SB_Filter.h>
#include <SB_Filter1D.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================

#include "MR_Obj.h"
float getNSigma(MRNoiseModel &self,int ind) {
return self.NSigma[ind];
}

void setNSigma(MRNoiseModel &self, int ind, float value) {
 self.NSigma[ind]=value;
}
float getTabEps(MRNoiseModel &self,int ind) {
return self.TabEps[ind];
}

void setTabEps(MRNoiseModel &self, int ind, float value) {
 self.TabEps[ind]=value;
}

boost::python::list band_to_scale(MultiResol &self, int ind) {
   int scale=0;
   details which_detail;
   boost::python::list retn;
   self.band_to_scale(ind,scale,which_detail); 
   retn.append(scale);
   retn.append(which_detail);
   return retn;
}




namespace  {

struct SubBand1D_Wrapper: SubBand1D
{
    SubBand1D_Wrapper(PyObject* py_self_, const SubBand1D& p0):
        SubBand1D(p0), py_self(py_self_) {}

    SubBand1D_Wrapper(PyObject* py_self_):
        SubBand1D(), py_self(py_self_) {}

    void transform(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "transform", p0, p1, p2, p3);
    }

    void default_transform(int p0, float* p1, float* p2, float* p3) {
        SubBand1D::transform(p0, p1, p2, p3);
    }

    void recons(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "recons", p0, p1, p2, p3);
    }

    void default_recons(int p0, float* p1, float* p2, float* p3) {
        SubBand1D::recons(p0, p1, p2, p3);
    }

    void noise_transform(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "noise_transform", p0, p1, p2, p3);
    }

    void default_noise_transform(int p0, float* p1, float* p2, float* p3) {
        SubBand1D::noise_transform(p0, p1, p2, p3);
    }

    void transform(int p0, float* p1, float* p2, float* p3, int p4) {
        call_method< void >(py_self, "transform", p0, p1, p2, p3, p4);
    }

    void default_transform(int p0, float* p1, float* p2, float* p3, int p4) {
        SubBand1D::transform(p0, p1, p2, p3, p4);
    }

    void recons(int p0, float* p1, float* p2, float* p3, int p4) {
        call_method< void >(py_self, "recons", p0, p1, p2, p3, p4);
    }

    void default_recons(int p0, float* p1, float* p2, float* p3, int p4) {
        SubBand1D::recons(p0, p1, p2, p3, p4);
    }

    PyObject* py_self;
};

struct SubBandFilter_Wrapper: SubBandFilter
{
    SubBandFilter_Wrapper(PyObject* py_self_, const SubBandFilter& p0):
        SubBandFilter(p0), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, type_sb_filter p0):
        SubBandFilter(p0), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, type_sb_filter p0, sb_type_norm p1):
        SubBandFilter(p0, p1), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, FilterAnaSynt& p0):
        SubBandFilter(p0), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, FilterAnaSynt& p0, sb_type_norm p1):
        SubBandFilter(p0, p1), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, char* p0):
        SubBandFilter(p0), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, char* p0, sb_type_norm p1):
        SubBandFilter(p0, p1), py_self(py_self_) {}

    SubBandFilter_Wrapper(PyObject* py_self_, float* p0, int p1, float* p2, int p3, sb_type_norm p4):
        SubBandFilter(p0, p1, p2, p3, p4), py_self(py_self_) {}

    void transform(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "transform", p0, p1, p2, p3);
    }

    void default_transform(int p0, float* p1, float* p2, float* p3) {
        SubBandFilter::transform(p0, p1, p2, p3);
    }

    void recons(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "recons", p0, p1, p2, p3);
    }

    void default_recons(int p0, float* p1, float* p2, float* p3) {
        SubBandFilter::recons(p0, p1, p2, p3);
    }

    void noise_transform(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "noise_transform", p0, p1, p2, p3);
    }

    void default_noise_transform(int p0, float* p1, float* p2, float* p3) {
        SubBandFilter::noise_transform(p0, p1, p2, p3);
    }

    void transform(int p0, float* p1, float* p2, float* p3, int p4) {
        call_method< void >(py_self, "transform", p0, p1, p2, p3, p4);
    }

    void default_transform(int p0, float* p1, float* p2, float* p3, int p4) {
        SubBandFilter::transform(p0, p1, p2, p3, p4);
    }

    void recons(int p0, float* p1, float* p2, float* p3, int p4) {
        call_method< void >(py_self, "recons", p0, p1, p2, p3, p4);
    }

    void default_recons(int p0, float* p1, float* p2, float* p3, int p4) {
        SubBandFilter::recons(p0, p1, p2, p3, p4);
    }

    PyObject* py_self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SubBandFilter_transform_overloads_2_3, transform, 2, 3)

struct UndecSubBandFilter_Wrapper: UndecSubBandFilter
{
    UndecSubBandFilter_Wrapper(PyObject* py_self_, const UndecSubBandFilter& p0):
        UndecSubBandFilter(p0), py_self(py_self_) {}

    UndecSubBandFilter_Wrapper(PyObject* py_self_):
        UndecSubBandFilter(), py_self(py_self_) {}

    UndecSubBandFilter_Wrapper(PyObject* py_self_, type_undec_filter p0):
        UndecSubBandFilter(p0), py_self(py_self_) {}

    void transform(int p0, float* p1, float* p2, float* p3, int p4) {
        call_method< void >(py_self, "transform", p0, p1, p2, p3, p4);
    }

    void default_transform(int p0, float* p1, float* p2, float* p3, int p4) {
        UndecSubBandFilter::transform(p0, p1, p2, p3, p4);
    }

    void recons(int p0, float* p1, float* p2, float* p3, int p4) {
        call_method< void >(py_self, "recons", p0, p1, p2, p3, p4);
    }

    void default_recons(int p0, float* p1, float* p2, float* p3, int p4) {
        UndecSubBandFilter::recons(p0, p1, p2, p3, p4);
    }

    void transform(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "transform", p0, p1, p2, p3);
    }

    void default_transform(int p0, float* p1, float* p2, float* p3) {
        SubBand1D::transform(p0, p1, p2, p3);
    }

    void recons(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "recons", p0, p1, p2, p3);
    }

    void default_recons(int p0, float* p1, float* p2, float* p3) {
        SubBand1D::recons(p0, p1, p2, p3);
    }

    void noise_transform(int p0, float* p1, float* p2, float* p3) {
        call_method< void >(py_self, "noise_transform", p0, p1, p2, p3);
    }

    void default_noise_transform(int p0, float* p1, float* p2, float* p3) {
        SubBand1D::noise_transform(p0, p1, p2, p3);
    }

    PyObject* py_self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SubBand2D_transform2d_overloads_1_6, transform2d, 1, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SubBand2D_recons2d_overloads_1_6, recons2d, 1, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_size_scale_nl_overloads_1_2, size_scale_nl, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_size_scale_nc_overloads_1_2, size_scale_nc, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_extract_scale_overloads_1_2, extract_scale, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_insert_scale_overloads_2_3, insert_scale, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_norm_overloads_0_2, norm, 0, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_alloc_overloads_4_8, alloc, 4, 8)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_computeNbBand_overloads_4_5, computeNbBand, 4, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_transform_overloads_2_3, transform, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_recons_overloads_1_2, recons, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_rec_adjoint_overloads_1_3, rec_adjoint, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_scale_norm_overloads_1_2, scale_norm, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_gauss_detect_level_overloads_1_3, gauss_detect_level, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_threshold_overloads_2_3, threshold, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultiResol_threshold_overloads_0_3, threshold, 0, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(ima_to_ortho_trans_overloads_2_3, ima_to_ortho_trans, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(ortho_trans_to_ima_overloads_2_3, ortho_trans_to_ima, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_alloc_overloads_5_9, alloc, 5, 9)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_prob_overloads_1_2, prob, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_prob_noise_overloads_1_2, prob_noise, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_kill_event_overloads_2_3, kill_event, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_threshold_overloads_1_2, threshold, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_weight_snr_overloads_1_2, weight_snr, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MRNoiseModel_weight_invsnr_overloads_1_2, weight_invsnr, 1, 2)


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(isap_great3_sparse2d)
{
    enum_< type_border >("type_border")
        .value("I_CONT", I_CONT)
        .value("I_ZERO", I_ZERO)
        .value("I_MIRROR", I_MIRROR)
        .value("I_PERIOD", I_PERIOD)
    ;

    enum_< type_filter >("type_filter")
        .value("FILTER_SOFT_THRESHOLD", FILTER_SOFT_THRESHOLD)
        .value("FILTER_MULTI_RES_WIENER", FILTER_MULTI_RES_WIENER)
        .value("FILTER_BAYES_BESSEL", FILTER_BAYES_BESSEL)
        .value("FILTER_MULTI_SOFT_MAD", FILTER_MULTI_SOFT_MAD)
        .value("FILTER_WAVELET_CONSTRAINT", FILTER_WAVELET_CONSTRAINT)
        .value("FILTER_ITER_ADJOINT_COEFF", FILTER_ITER_ADJOINT_COEFF)
        .value("FILTER_MULTI_HARD_MAD", FILTER_MULTI_HARD_MAD)
        .value("FILTER_THRESHOLD", FILTER_THRESHOLD)
        .value("FILTER_HIERARCHICAL_TRESHOLD", FILTER_HIERARCHICAL_TRESHOLD)
        .value("FILTER_BIVARIATE_SHRINKAGE", FILTER_BIVARIATE_SHRINKAGE)
        .value("FILTER_TV_CONSTRAINT", FILTER_TV_CONSTRAINT)
        .value("FILTER_ITER_THRESHOLD", FILTER_ITER_THRESHOLD)
        .value("FILTER_HIERARCHICAL_WIENER", FILTER_HIERARCHICAL_WIENER)
    ;

    enum_< sb_type_norm >("sb_type_norm")
        .value("NORM_L2", NORM_L2)
        .value("NORM_L1", NORM_L1)
    ;

    enum_< type_sb_filter >("type_sb_filter")
        .value("SB_UNKNOWN", SB_UNKNOWN)
        .value("F_4_4", F_4_4)
        .value("F_BI2HAAR", F_BI2HAAR)
        .value("F_MALLAT_9_7", F_MALLAT_9_7)
        .value("F_LEMARIE_1", F_LEMARIE_1)
        .value("F_USER", F_USER)
        .value("F_LEMARIE_3", F_LEMARIE_3)
        .value("F_5_3", F_5_3)
        .value("F_LEMARIE_5", F_LEMARIE_5)
        .value("F_5_3_DIV", F_5_3_DIV)
        .value("F_MALLAT_7_9", F_MALLAT_7_9)
        .value("F_3_5", F_3_5)
        .value("F_HAAR", F_HAAR)
        .value("F_ODEGARD_7_9", F_ODEGARD_7_9)
        .value("F_DAUBE_4", F_DAUBE_4)
        .value("F_BI4HAAR", F_BI4HAAR)
    ;

    class_< FilterAnaSynt >("FilterAnaSynt", init<  >())
        .def(init< const FilterAnaSynt& >())
        .def(init< type_sb_filter >())
        .def(init< char* >())
        .def_readwrite("Verbose", &FilterAnaSynt::Verbose)
        .def_readwrite("FilterFileName", &FilterAnaSynt::FilterFileName)
        .def("type_filter", &FilterAnaSynt::type_filter)
        .def("norm", &FilterAnaSynt::norm)
        .def("size_analysis", &FilterAnaSynt::size_analysis)
        .def("size_synthesis", &FilterAnaSynt::size_synthesis)
        .def("start_analysis", &FilterAnaSynt::start_analysis)
        .def("start_synthesis", &FilterAnaSynt::start_synthesis)
        .def("alloc", &FilterAnaSynt::alloc)
    ;

    def("StringSBFilter", &StringSBFilter);
    enum_< type_undec_filter >("type_undec_filter")
        .value("U_B3SPLINE", U_B3SPLINE)
        .value("U_HAAR_B3S_POS", U_HAAR_B3S_POS)
        .value("U_HAAR_B3S", U_HAAR_B3S)
        .value("U_B3SPLINE_2", U_B3SPLINE_2)
        .value("U_B2SPLINE", U_B2SPLINE)
    ;

    enum_< type_lift >("type_lift")
        .value("TL_UNKNOWN", TL_UNKNOWN)
        .value("TL_INT_HAAR", TL_INT_HAAR)
        .value("TL_INT_CDF", TL_INT_CDF)
        .value("TL_CDF", TL_CDF)
        .value("TL_MEDIAN", TL_MEDIAN)
        .value("TL_INT_4_2", TL_INT_4_2)
        .value("TL_INT_F79", TL_INT_F79)
        .value("TL_F79", TL_F79)
    ;

    class_< SubBand1D, SubBand1D_Wrapper >("SubBand1D", init<  >())
        .def(init< const SubBand1D& >())
        .def_readwrite("SubSample_H_Even", &SubBand1D::SubSample_H_Even)
        .def_readwrite("SubSample_G_Odd", &SubBand1D::SubSample_G_Odd)
        .def_readwrite("DistPix", &SubBand1D::DistPix)
        .def_readwrite("Border", &SubBand1D::Border)
        .def("transform", (void (SubBand1D::*)(int, float*, float*, float*) )&SubBand1D::transform, (void (SubBand1D_Wrapper::*)(int, float*, float*, float*))&SubBand1D_Wrapper::default_transform)
        .def("recons", (void (SubBand1D::*)(int, float*, float*, float*) )&SubBand1D::recons, (void (SubBand1D_Wrapper::*)(int, float*, float*, float*))&SubBand1D_Wrapper::default_recons)
        .def("noise_transform", &SubBand1D::noise_transform, &SubBand1D_Wrapper::default_noise_transform)
        .def("transform", (void (SubBand1D::*)(int, float*, float*, float*, int) )&SubBand1D::transform, (void (SubBand1D_Wrapper::*)(int, float*, float*, float*, int))&SubBand1D_Wrapper::default_transform)
        .def("recons", (void (SubBand1D::*)(int, float*, float*, float*, int) )&SubBand1D::recons, (void (SubBand1D_Wrapper::*)(int, float*, float*, float*, int))&SubBand1D_Wrapper::default_recons)
        .def("test_index", &SubBand1D::test_index)
    ;

    class_< SubBandFilter, bases< SubBand1D > , SubBandFilter_Wrapper >("SubBandFilter", init< const SubBandFilter& >())
        .def(init< type_sb_filter >())
        .def(init< type_sb_filter, sb_type_norm >())
        .def(init< FilterAnaSynt& >())
        .def(init< FilterAnaSynt&, sb_type_norm >())
        .def(init< char* >())
        .def(init< char*, sb_type_norm >())
        .def(init< float*, int, float*, int, sb_type_norm >())
        .def("transform", (void (SubBandFilter::*)(int, float*, float*, float*) )&SubBandFilter::transform, (void (SubBandFilter_Wrapper::*)(int, float*, float*, float*))&SubBandFilter_Wrapper::default_transform)
        .def("recons", (void (SubBandFilter::*)(int, float*, float*, float*) )&SubBandFilter::recons, (void (SubBandFilter_Wrapper::*)(int, float*, float*, float*))&SubBandFilter_Wrapper::default_recons)
        .def("noise_transform", (void (SubBandFilter::*)(int, float*, float*, float*) )&SubBandFilter::noise_transform, (void (SubBandFilter_Wrapper::*)(int, float*, float*, float*))&SubBandFilter_Wrapper::default_noise_transform)
        .def("transform", (void (SubBandFilter::*)(int, float*, float*, float*, int) )&SubBandFilter::transform, (void (SubBandFilter_Wrapper::*)(int, float*, float*, float*, int))&SubBandFilter_Wrapper::default_transform)
        .def("recons", (void (SubBandFilter::*)(int, float*, float*, float*, int) )&SubBandFilter::recons, (void (SubBandFilter_Wrapper::*)(int, float*, float*, float*, int))&SubBandFilter_Wrapper::default_recons)
        .def("convol_h0", &SubBandFilter::convol_h0)
        .def("convol_g0", &SubBandFilter::convol_g0)
        .def("noise_convol_h0", &SubBandFilter::noise_convol_h0)
        .def("noise_convol_g0", &SubBandFilter::noise_convol_g0)
        .def("convol_h1", &SubBandFilter::convol_h1)
        .def("convol_g1", &SubBandFilter::convol_g1)
        .def("convol_filter", &SubBandFilter::convol_filter)
        .def("rec_convol_filter", &SubBandFilter::rec_convol_filter)
        .def("transform", (void (SubBandFilter::*)(to_array<float,true>&, to_array<float,true>&, to_array<float,true>*) )&SubBandFilter::transform, SubBandFilter_transform_overloads_2_3())
        .def("recons", (void (SubBandFilter::*)(to_array<float,true>&, to_array<float,true>&, to_array<float,true>&) )&SubBandFilter::recons)
    ;

    class_< UndecSubBandFilter, bases< SubBand1D > , UndecSubBandFilter_Wrapper >("UndecSubBandFilter", init<  >())
        .def(init< const UndecSubBandFilter& >())
        .def(init< type_undec_filter >())
        .def_readwrite("Verbose", &UndecSubBandFilter::Verbose)
        .def("transform", (void (UndecSubBandFilter::*)(int, float*, float*, float*, int) )&UndecSubBandFilter::transform, (void (UndecSubBandFilter_Wrapper::*)(int, float*, float*, float*, int))&UndecSubBandFilter_Wrapper::default_transform)
        .def("recons", (void (UndecSubBandFilter::*)(int, float*, float*, float*, int) )&UndecSubBandFilter::recons, (void (UndecSubBandFilter_Wrapper::*)(int, float*, float*, float*, int))&UndecSubBandFilter_Wrapper::default_recons)
        .def("transform", (void (SubBand1D::*)(int, float*, float*, float*) )&SubBand1D::transform, (void (UndecSubBandFilter_Wrapper::*)(int, float*, float*, float*))&UndecSubBandFilter_Wrapper::default_transform)
        .def("recons", (void (SubBand1D::*)(int, float*, float*, float*) )&SubBand1D::recons, (void (UndecSubBandFilter_Wrapper::*)(int, float*, float*, float*))&UndecSubBandFilter_Wrapper::default_recons)
        .def("noise_transform", &SubBand1D::noise_transform, &UndecSubBandFilter_Wrapper::default_noise_transform)
    ;

    def("StringLSTransform", &StringLSTransform);
    def("StringUndecFilter", &StringUndecFilter);
    class_< SubBand2D >("SubBand2D", init< const SubBand2D& >())
        .def(init< SubBand1D& >())
        .def(init< SubBand1D&, SubBand1D& >())
        .def_readwrite("Line_SubSample_H_Even", &SubBand2D::Line_SubSample_H_Even)
        .def_readwrite("Line_SubSample_G_Odd", &SubBand2D::Line_SubSample_G_Odd)
        .def_readwrite("Col_SubSample_H_Even", &SubBand2D::Col_SubSample_H_Even)
        .def_readwrite("Col_SubSample_G_Odd", &SubBand2D::Col_SubSample_G_Odd)
        .def("get_subband_method", &SubBand2D::get_subband_method, return_internal_reference< 1 >())
        .def("get_subband_method_line", &SubBand2D::get_subband_method_line, return_internal_reference< 1 >())
        .def("get_subband_method_col", &SubBand2D::get_subband_method_col, return_internal_reference< 1 >())
        .def("transform2d", (void (SubBand2D::*)(to_array<float,false>&, Bool, to_array<float,false>*, to_array<float,false>*, to_array<float,false>*, to_array<float,false>*) )&SubBand2D::transform2d, SubBand2D_transform2d_overloads_1_6())
        .def("recons2d", (void (SubBand2D::*)(to_array<float,false>&, Bool, to_array<float,false>*, to_array<float,false>*, to_array<float,false>*, to_array<float,false>*) )&SubBand2D::recons2d, SubBand2D_recons2d_overloads_1_6())
        .def("transform2d", (void (SubBand2D::*)(to_array<float,false>&, to_array<float,false>&, to_array<float,false>&, to_array<float,false>&, to_array<float,false>&, int) )&SubBand2D::transform2d)
        .def("recons2d", (void (SubBand2D::*)(to_array<float,false>&, to_array<float,false>&, to_array<float,false>&, to_array<float,false>&, to_array<float,false>&, int) )&SubBand2D::recons2d)
    ;

    enum_< details >("details")
        .value("I_SMOOTH", I_SMOOTH)
        .value("D_DIAGONAL", D_DIAGONAL)
        .value("D_HALF_RESOL", D_HALF_RESOL)
        .value("D_HORIZONTAL", D_HORIZONTAL)
        .value("D_RESOL", D_RESOL)
        .value("D_NULL", D_NULL)
        .value("D_VERTICAL", D_VERTICAL)
    ;

    enum_< type_transform >("type_transform")
        .value("TM_TO_SEMI_PYR", TM_TO_SEMI_PYR)
        .value("TM_PAVE_MEDIAN", TM_PAVE_MEDIAN)
        .value("TO_SEMI_PYR", TO_SEMI_PYR)
        .value("TO_HAAR", TO_HAAR)
        .value("TM_PYR_SCALING_FUNCTION", TM_PYR_SCALING_FUNCTION)
        .value("TO_PYR_BSPLINE", TO_PYR_BSPLINE)
        .value("TM_PYR_MINMAX", TM_PYR_MINMAX)
        .value("TO_FEAUVEAU", TO_FEAUVEAU)
        .value("T_UNDEFINED", T_UNDEFINED)
        .value("TO_PAVE_FEAUVEAU", TO_PAVE_FEAUVEAU)
        .value("TO_UNDECIMATED_MALLAT", TO_UNDECIMATED_MALLAT)
        .value("TO_PYR_LINEAR", TO_PYR_LINEAR)
        .value("TM_MIN_MAX", TM_MIN_MAX)
        .value("TO_PAVE_FFT", TO_PAVE_FFT)
        .value("TM_PAVE_MINMAX", TM_PAVE_MINMAX)
        .value("TM_PYR_MEDIAN", TM_PYR_MEDIAN)
        .value("TO_DIADIC_HAAR", TO_DIADIC_HAAR)
        .value("TO_PAVE_LINEAR", TO_PAVE_LINEAR)
        .value("TO_LIFTING", TO_LIFTING)
        .value("TM_TO_PYR", TM_TO_PYR)
        .value("TO_DIV_2", TO_DIV_2)
        .value("TO_DIV_1", TO_DIV_1)
        .value("TO_UNDECIMATED_NON_ORTHO", TO_UNDECIMATED_NON_ORTHO)
        .value("TO_PYR_FFT_DIFF_RESOL", TO_PYR_FFT_DIFF_RESOL)
        .value("TO_LC", TO_LC)
        .value("TO_PYR_MEYER", TO_PYR_MEYER)
        .value("TO_PAVE_HAAR", TO_PAVE_HAAR)
        .value("TO_PYR_FFT_DIFF_SQUARE", TO_PYR_FFT_DIFF_SQUARE)
        .value("TM_PYR_LAPLACIAN", TM_PYR_LAPLACIAN)
        .value("TO_PYR_MEYER_ISOTROP", TO_PYR_MEYER_ISOTROP)
        .value("TO_PAVE_BSPLINE", TO_PAVE_BSPLINE)
        .value("TO_MALLAT", TO_MALLAT)
        .value("TC_FCT", TC_FCT)
        .value("TO_DIADIC_MALLAT", TO_DIADIC_MALLAT)
    ;

    enum_< set_transform >("set_transform")
        .value("TRANSF_SEMIPYR", TRANSF_SEMIPYR)
        .value("TRANSF_PYR", TRANSF_PYR)
        .value("TRANSF_UNDECIMATED_MALLAT", TRANSF_UNDECIMATED_MALLAT)
        .value("TRANSF_PAVE", TRANSF_PAVE)
        .value("TRANSF_MALLAT", TRANSF_MALLAT)
        .value("TRANSF_FEAUVEAU", TRANSF_FEAUVEAU)
        .value("TRANSF_DIADIC_MALLAT", TRANSF_DIADIC_MALLAT)
        .value("S_UNDEFINED", S_UNDEFINED)
    ;

    def("StringTransform", &StringTransform);
    enum_< type_noise >("type_noise")
        .value("NOISE_EVENT_POISSON", NOISE_EVENT_POISSON)
        .value("NOISE_CORREL", NOISE_CORREL)
        .value("NOISE_GAUSSIAN", NOISE_GAUSSIAN)
        .value("NOISE_SPECKLE", NOISE_SPECKLE)
        .value("NOISE_NON_UNI_ADD", NOISE_NON_UNI_ADD)
        .value("NOISE_UNDEFINED", NOISE_UNDEFINED)
        .value("NOISE_UNI_UNDEFINED", NOISE_UNI_UNDEFINED)
        .value("NOISE_POISSON", NOISE_POISSON)
        .value("NOISE_NON_UNI_MULT", NOISE_NON_UNI_MULT)
        .value("NOISE_MULTI", NOISE_MULTI)
        .value("NOISE_GAUSS_POISSON", NOISE_GAUSS_POISSON)
    ;

    class_< MultiResol >("MultiResol", init<  >())
        .def(init< const MultiResol& >())
        .def(init< int, int, int, type_transform, char* >())
        .def_readwrite("FCT_Nbr_Angle", &MultiResol::FCT_Nbr_Angle)
        .def_readwrite("ModifiedATWT", &MultiResol::ModifiedATWT)
        .def_readwrite("Verbose", &MultiResol::Verbose)
        .def_readwrite("Type_Transform", &MultiResol::Type_Transform)
        .def_readwrite("Set_Transform", &MultiResol::Set_Transform)
        .def_readwrite("FormatInputImag", &MultiResol::FormatInputImag)
        .def_readwrite("Border", &MultiResol::Border)
        .def_readwrite("MedianWindowSize", &MultiResol::MedianWindowSize)
        .def_readwrite("Fc", &MultiResol::Fc)
        .def_readwrite("Nbr_Iter", &MultiResol::Nbr_Iter)
        .def_readwrite("ExactPyrRec", &MultiResol::ExactPyrRec)
        .def_readwrite("SigmaNoise", &MultiResol::SigmaNoise)
        .def_readwrite("TMTO_NSigBSpline", &MultiResol::TMTO_NSigBSpline)
        .def_readwrite("TMTO_NSigMedian", &MultiResol::TMTO_NSigMedian)
        .def_readwrite("LiftingTrans", &MultiResol::LiftingTrans)
        .def_readwrite("SBFilter", &MultiResol::SBFilter)
        .def_readwrite("U_Filter", &MultiResol::U_Filter)
        .def_readwrite("TypeNorm", &MultiResol::TypeNorm)
        .def_readwrite("NormHaar", &MultiResol::NormHaar)
        .def_readwrite("EdgeLineTransform", &MultiResol::EdgeLineTransform)
        .def_readwrite("LC", &MultiResol::LC)
        .def_readwrite("LC_NbrScaleLine", &MultiResol::LC_NbrScaleLine)
        .def_readwrite("LC_NbrScaleCol", &MultiResol::LC_NbrScaleCol)
        .def("filter_bank", &MultiResol::filter_bank, return_internal_reference< 1 >())
        .def("filter_bank_line", &MultiResol::filter_bank_line, return_internal_reference< 1 >())
        .def("filter_bank_column", &MultiResol::filter_bank_column, return_internal_reference< 1 >())
        .def("undec_filter_bank", &MultiResol::undec_filter_bank, return_internal_reference< 1 >())
        .def("nbr_undec_scale", &MultiResol::nbr_undec_scale)
        .def("nbr_scale", &MultiResol::nbr_scale)
        .def("nbr_band", &MultiResol::nbr_band)
        .def("nbr_band_per_resol", (int (MultiResol::*)() const)&MultiResol::nbr_band_per_resol)
        .def("nbr_band_per_resol", (int (MultiResol::*)(int) const)&MultiResol::nbr_band_per_resol)
        .def("nbr_mr_coeff", &MultiResol::nbr_mr_coeff)
        .def("scale_to_band", &MultiResol::scale_to_band)
        .def("stb", &MultiResol::stb)
        .def("nbr_coeff_in_band", &MultiResol::nbr_coeff_in_band)
        .def("size_band_nl", &MultiResol::size_band_nl)
        .def("size_band_nc", &MultiResol::size_band_nc)
        .def("size_scale_nl", &MultiResol::size_scale_nl, MultiResol_size_scale_nl_overloads_1_2())
        .def("size_scale_nc", &MultiResol::size_scale_nc, MultiResol_size_scale_nc_overloads_1_2())
        .def("size_ima_nl", &MultiResol::size_ima_nl)
        .def("size_ima_nc", &MultiResol::size_ima_nc)
        .def("pos_coeff", &MultiResol::pos_coeff)
        .def("extract_scale", &MultiResol::extract_scale, MultiResol_extract_scale_overloads_1_2())
        .def("extract_band", &MultiResol::extract_band)
        .def("insert_scale", &MultiResol::insert_scale, MultiResol_insert_scale_overloads_2_3())
        .def("insert_band", &MultiResol::insert_band)
        .def("norm", &MultiResol::norm, MultiResol_norm_overloads_0_2())
        .def("alloc", (void (MultiResol::*)(int, int, int, type_transform, char*) )&MultiResol::alloc)
        .def("alloc", (void (MultiResol::*)(int, int, int, type_transform, FilterAnaSynt*, sb_type_norm, int, type_undec_filter) )&MultiResol::alloc, MultiResol_alloc_overloads_4_8())
        .def("computeNbBand", &MultiResol::computeNbBand, MultiResol_computeNbBand_overloads_4_5())
        .def("free", &MultiResol::free)
        .def("band", &MultiResol::band, return_internal_reference< 1 >())
        .def("tabband", &MultiResol::tabband, return_internal_reference< 1 >())
        .def("transform", (void (MultiResol::*)(to_array<float,false>&, type_border, Bool) )&MultiResol::transform, MultiResol_transform_overloads_2_3())
        .def("compute_mod_phase", &MultiResol::compute_mod_phase)
        .def("transform", (void (MultiResol::*)(to_array<float,false>&) )&MultiResol::transform)
        .def("recons", &MultiResol::recons, MultiResol_recons_overloads_1_2())
        .def("rec_adjoint", &MultiResol::rec_adjoint, MultiResol_rec_adjoint_overloads_1_3())
        .def("scale_norm", &MultiResol::scale_norm, MultiResol_scale_norm_overloads_1_2())
        .def("band_norm", &MultiResol::band_norm)
        .def("read", (void (MultiResol::*)(char*) )&MultiResol::read)
        .def("write", (void (MultiResol::*)(char*) )&MultiResol::write)
        .def("read", (void (MultiResol::*)(char*, int) )&MultiResol::read)
        .def("write", (void (MultiResol::*)(char*, int) )&MultiResol::write)
        .def("read_band", &MultiResol::read_band)
        .def("write_band", &MultiResol::write_band)
        .def("gauss_detect_level", &MultiResol::gauss_detect_level, MultiResol_gauss_detect_level_overloads_1_3())
        .def("threshold", (void (MultiResol::*)(int, float, Bool) )&MultiResol::threshold, MultiResol_threshold_overloads_2_3())
        .def("threshold", (void (MultiResol::*)(float, float, Bool) )&MultiResol::threshold, MultiResol_threshold_overloads_0_3())
        .def("print_info", &MultiResol::print_info)
        .def("band_to_scale",&band_to_scale)
    ;

    def("ima_to_ortho_trans", &ima_to_ortho_trans, ima_to_ortho_trans_overloads_2_3());
    def("ortho_trans_to_ima", &ortho_trans_to_ima, ortho_trans_to_ima_overloads_2_3());
    enum_< type_sigma_method >("type_sigma_method")
        .value("SIGMA_CLIPIMA", SIGMA_CLIPIMA)
        .value("SIGMA_SUPPORT", SIGMA_SUPPORT)
        .value("SIGMA_MEDIAN", SIGMA_MEDIAN)
        .value("SIGMA_BSPLINE", SIGMA_BSPLINE)
        .value("SIGM_UNDEFINED", SIGM_UNDEFINED)
    ;

    class_< MRNoiseModel >("MRNoiseModel", init<  >())
        .def(init< const MRNoiseModel& >())
        .def(init< type_noise, int, int, int, type_transform >())
        .def(init< type_noise, MultiResol& >())
        .def_readwrite("TypeThreshold", &MRNoiseModel::TypeThreshold)
        .def_readwrite("U_Filter", &MRNoiseModel::U_Filter)
        .def_readwrite("TypeNorm", &MRNoiseModel::TypeNorm)
        .def_readwrite("TransImag", &MRNoiseModel::TransImag)
        .def_readwrite("NiterSigmaClip", &MRNoiseModel::NiterSigmaClip)
        .def_readwrite("SizeBlockSigmaNoise", &MRNoiseModel::SizeBlockSigmaNoise)
        .def_readwrite("NoCompDistrib", &MRNoiseModel::NoCompDistrib)
        .def_readwrite("EventImage", &MRNoiseModel::EventImage)
        .def_readwrite("DilateSupport", &MRNoiseModel::DilateSupport)
        .def_readwrite("GetEgde", &MRNoiseModel::GetEgde)
        .def_readwrite("SupIsol", &MRNoiseModel::SupIsol)
        .def_readwrite("MinEvent", &MRNoiseModel::MinEvent)
        .def_readwrite("OnlyPositivDetect", &MRNoiseModel::OnlyPositivDetect)
        .def_readwrite("CFewEventPoisson2d", &MRNoiseModel::CFewEventPoisson2d)
        .def_readwrite("CFewEvent2d", &MRNoiseModel::CFewEvent2d)
        .def_readwrite("MinEventNumber", &MRNoiseModel::MinEventNumber)
        .def_readwrite("Event_Image", &MRNoiseModel::Event_Image)
        .def_readwrite("UseRmsMap", &MRNoiseModel::UseRmsMap)
        .def_readwrite("RmsMap", &MRNoiseModel::RmsMap)
        .def_readwrite("GetRmsCoeffGauss", &MRNoiseModel::GetRmsCoeffGauss)
        .def_readwrite("MR_Data_Event", &MRNoiseModel::MR_Data_Event)
        .def_readwrite("BadPixel", &MRNoiseModel::BadPixel)
        .def_readwrite("BadPixelVal", &MRNoiseModel::BadPixelVal)
        .def_readwrite("FirstDectectScale", &MRNoiseModel::FirstDectectScale)
        .def_readwrite("SigmaDetectionMethod", &MRNoiseModel::SigmaDetectionMethod)
        .def_readwrite("CCD_Gain", &MRNoiseModel::CCD_Gain)
        .def_readwrite("CCD_ReadOutSigma", &MRNoiseModel::CCD_ReadOutSigma)
        .def_readwrite("CCD_ReadOutMean", &MRNoiseModel::CCD_ReadOutMean)
        .def_readwrite("SigmaNoise", &MRNoiseModel::SigmaNoise)
        .def_readwrite("CorrelNoiseMap", &MRNoiseModel::CorrelNoiseMap)
        .def_readwrite("CSpeckle", &MRNoiseModel::CSpeckle)
        .def_readwrite("SigmaApprox", &MRNoiseModel::SigmaApprox)
        .def_readwrite("GradientAnalysis", &MRNoiseModel::GradientAnalysis)
        .def_readwrite("PoissonFisz", &MRNoiseModel::PoissonFisz)
        .def_readwrite("MadEstimFromCenter", &MRNoiseModel::MadEstimFromCenter)
        .def("filter_bank", &MRNoiseModel::filter_bank, return_internal_reference< 1 >())
        .def("nbr_scale", &MRNoiseModel::nbr_scale)
        .def("nbr_band", &MRNoiseModel::nbr_band)
        .def("nbr_undec_scale", &MRNoiseModel::nbr_undec_scale)
        .def("nbr_coeff_in_band", &MRNoiseModel::nbr_coeff_in_band)
        .def("nl", &MRNoiseModel::nl)
        .def("nc", &MRNoiseModel::nc)
        .def("which_noise", &MRNoiseModel::which_noise)
        .def("type_trans", &MRNoiseModel::type_trans)
        .def("alloc", &MRNoiseModel::alloc, MRNoiseModel_alloc_overloads_5_9())
        .def("free", &MRNoiseModel::free)
        .def("model", (void (MRNoiseModel::*)(to_array<float,false>&) )&MRNoiseModel::model)
        .def("model", (void (MRNoiseModel::*)(to_array<float,false>&, MultiResol&) )&MRNoiseModel::model)
        .def("write_support_mr", &MRNoiseModel::write_support_mr)
        .def("write_support_ima", &MRNoiseModel::write_support_ima)
        .def("mr_obj", &MRNoiseModel::mr_obj)
        .def("write_in_few_event", &MRNoiseModel::write_in_few_event)
        .def("nsigma", &MRNoiseModel::nsigma)
        .def("signif", (Bool (MRNoiseModel::*)(float, int, int, int, float, float) )&MRNoiseModel::signif)
        .def("signif", (Bool (MRNoiseModel::*)(float, int, int, int) )&MRNoiseModel::signif)
        .def("signif", (Bool (MRNoiseModel::*)(float, int, int, int, to_array<float,true>&) )&MRNoiseModel::signif)
        .def("signif", (Bool (MRNoiseModel::*)(float, int, int, int, details) )&MRNoiseModel::signif)
        .def("prob", (float (MRNoiseModel::*)(float, int, int, int, details) )&MRNoiseModel::prob)
        .def("prob", (float (MRNoiseModel::*)(float, int, int, int) )&MRNoiseModel::prob)
        .def("prob", (void (MRNoiseModel::*)(MultiResol&, Bool) )&MRNoiseModel::prob, MRNoiseModel_prob_overloads_1_2())
        .def("prob_noise", (void (MRNoiseModel::*)(MultiResol&, Bool) )&MRNoiseModel::prob_noise, MRNoiseModel_prob_noise_overloads_1_2())
        .def("prob_noise", (float (MRNoiseModel::*)(float, int, int, int) )&MRNoiseModel::prob_noise)
        .def("prob_noise", (float (MRNoiseModel::*)(float, int, int, int, details) )&MRNoiseModel::prob_noise)
        .def("prob_signal", (float (MRNoiseModel::*)(float, int, int, int) )&MRNoiseModel::prob_signal)
        .def("prob_signal", (float (MRNoiseModel::*)(float, int, int, int, details) )&MRNoiseModel::prob_signal)
        .def("im_transform", &MRNoiseModel::im_transform)
        .def("im_invtransform", &MRNoiseModel::im_invtransform)
        .def("val_transform", &MRNoiseModel::val_transform)
        .def("val_invtransform", &MRNoiseModel::val_invtransform)
        .def("kill_isol", &MRNoiseModel::kill_isol)
        .def("dilate_support", &MRNoiseModel::dilate_support)
        .def("kill_event", &MRNoiseModel::kill_event, MRNoiseModel_kill_event_overloads_2_3())
        .def("set_support", &MRNoiseModel::set_support)
        .def("mod_support", &MRNoiseModel::mod_support)
        .def("refresh_support", &MRNoiseModel::refresh_support)
        .def("set_sigma", &MRNoiseModel::set_sigma)
        .def("threshold", &MRNoiseModel::threshold, MRNoiseModel_threshold_overloads_1_2())
        .def("weight_snr", &MRNoiseModel::weight_snr, MRNoiseModel_weight_snr_overloads_1_2())
        .def("weight_invsnr", &MRNoiseModel::weight_invsnr, MRNoiseModel_weight_invsnr_overloads_1_2())
        .def("prob_signal_few_event", &MRNoiseModel::prob_signal_few_event)
        .def("hierarchical_dilate_support", &MRNoiseModel::hierarchical_dilate_support)
        .def("set_old_poisson", &MRNoiseModel::set_old_poisson)
        .def("write_threshold", &MRNoiseModel::write_threshold)
        .def("old_poisson", &MRNoiseModel::old_poisson)
        .def("trace", &MRNoiseModel::trace)
        .def("__getNSigma",&getNSigma)
        .def("__setNSigma",&setNSigma)
        .def("__getTabEps",&getTabEps)
        .def("__setTabEps",&setTabEps)
    ;

}

