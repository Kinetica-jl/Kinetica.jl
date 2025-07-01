const PYTYPES = setdiff(getindex.(getfield.(getfield.(methods(ispy), :sig), :types), 2), [Any])
const PYCONVERT_NATIVE_TYPES_CACHE = Dict{PythonCall.C.PyPtr, Vector{Type}}()

function pyconvert_get_native_ruletypes(x)
    tptr = PythonCall.C.Py_Type(PythonCall.Core.getptr(x))
    get!(PYCONVERT_NATIVE_TYPES_CACHE , tptr) do
        pytype = PythonCall.pynew(PythonCall.Compat.incref(tptr))
        ans = setdiff(getfield.(PythonCall.Convert._pyconvert_get_rules(pytype), :type), PYTYPES)
        PythonCall.Core.pydel!(pytype)
        ans
    end
end

function pyconvert_native(x)
    ispy(x) || return x
    
    PythonCall.pyisinstance(x, pybuiltins.dict) && return pyconvert_dict(x)

    types = pyconvert_get_native_ruletypes(x)
    for T in types 
        ans = PythonCall.Convert.pytryconvert(T, x)
        PythonCall.Convert.pyconvert_isunconverted(ans) || return PythonCall.Convert.pyconvert_result(T, ans)
    end
    PythonCall.Convert.pyconvert_unconverted()
end

function pyconvert_dict(::Type{T}, d) where T
    jd = pyconvert(T, d)
    kk = pyconvert_native.(keys(jd))
    vv = [ispy(v) ? PythonCall.pyisinstance(v, pybuiltins.dict) ? pyconvert_dict(T, v) : pyconvert_native(v) : v for v in values(jd)]

    T(zip(kk, vv))
end

pyconvert_dict(d) = pyconvert_dict(Dict, d)