cdef extern from "pyrex_helper.hpp":
    void __helper_set_double "pyrex_helper::ref_set<double>" (double, double)
    void __helper_set_float "pyrex_helper::ref_set<float>" (float, float)
    void __helper_set_long "pyrex_helper::ref_set<int>" (int, int)
    void __helper_set_int "pyrex_helper::ref_set<int>" (int, int)
    void __helper_set_short "pyrex_helper::ref_set<short>" (short, short)
    void __helper_set_char "pyrex_helper::ref_set<char>" (char, char)
    object __helper_pystr_from_repr "pyrex_helper::pystr_from_repr" (void *)
