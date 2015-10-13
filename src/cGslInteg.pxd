
cdef extern from "gsl/gsl_math.h":
    
    # Definition of an arbitrary function with parameters
    ctypedef struct gsl_function:
        double (* function) (double x, void * params)
        void * params
    

cdef extern from "gsl/gsl_integration.h":
    
    
    #"""Re-define all functions, structures to be used in the python code. Mostly can copy these
    #   straight from the .h file"""
       
    # opaque handle: cython does not need to know the contents of this struct
    ctypedef struct gsl_integration_workspace
        
    gsl_integration_workspace * gsl_integration_workspace_alloc (const size_t n)
    
    void gsl_integration_workspace_free(gsl_integration_workspace* w)
    
    # qags = quadrature + adaptive + general (user-defined) integrand + singularities taken care of
    int gsl_integration_qags (const gsl_function * f,
                          double a, double b,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);    


