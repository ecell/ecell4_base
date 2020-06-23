#include <ecell4/core/PolygonHDF5Writer.hpp>
#include <ecell4/core/Polygon.hpp>

namespace ecell4
{

struct PolygonHDF5Traits
{
    struct h5_triangle_struct {
        double p1x;
        double p1y;
        double p1z;
        double p2x;
        double p2y;
        double p2z;
        double p3x;
        double p3y;
        double p3z;
    };

    static H5::CompType get_triangle_comp_type()
    {
        H5::CompType h5_triangle_comp_type(sizeof(h5_triangle_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(h5_triangle_comp_type.getId(), #member,\
                HOFFSET(h5_triangle_struct, member), type.getId())
        INSERT_MEMBER(p1x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p1y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p1z, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p2x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p2y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p2z, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p3x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p3y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p3z, H5::PredType::NATIVE_DOUBLE);
#undef INSERT_MEMBER
        return h5_triangle_comp_type;
    }
};

void save_polygon_hdf5(const Polygon& p, H5::Group* root)
{
    using traits_type        = PolygonHDF5Traits;
    using h5_triangle_struct = traits_type::h5_triangle_struct;

    // ----------------------------------------------------------------------
    // convert face data into the simpler type, h5_triangle_structs

    const auto triangles             = p.triangles();
    const unsigned int num_triangles = triangles.size();

    std::unique_ptr<h5_triangle_struct[]>
        h5_triangle_table(new h5_triangle_struct[num_triangles]);
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
        const auto triangle = triangles.at(i);
        h5_triangle_table[i].p1x = triangle.vertex_at(0)[0];
        h5_triangle_table[i].p1y = triangle.vertex_at(0)[1];
        h5_triangle_table[i].p1z = triangle.vertex_at(0)[2];
        h5_triangle_table[i].p2x = triangle.vertex_at(1)[0];
        h5_triangle_table[i].p2y = triangle.vertex_at(1)[1];
        h5_triangle_table[i].p2z = triangle.vertex_at(1)[2];
        h5_triangle_table[i].p3x = triangle.vertex_at(2)[0];
        h5_triangle_table[i].p3y = triangle.vertex_at(2)[1];
        h5_triangle_table[i].p3z = triangle.vertex_at(2)[2];
    }

    // ----------------------------------------------------------------------
    // save triangles

    const int     RANK   = 1;
    const hsize_t dim1[] = {num_triangles};
    H5::DataSpace dataspace(RANK, dim1);
    std::unique_ptr<H5::DataSet> dataset(new H5::DataSet(root->createDataSet(
            "triangles", traits_type::get_triangle_comp_type(), dataspace)));
    dataset->write(h5_triangle_table.get(), dataset->getDataType());

    // ----------------------------------------------------------------------
    // save edge lengths
    const Real3 edge_lengths = p.edge_lengths();
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);

    // ----------------------------------------------------------------------
    // Since `Polygon` does not allow to remove/insert new triangle, the
    // triangles saved here should have been assigned all at once.
    // To make Identifier consistent, we need to call `polygon.assign` while
    // loading. If we did that, we don't need to save the corresponding IDs
    // because polygon->assign() function assigns the same IDs, in the same order.

    return;
}

void load_polygon_hdf5(const H5::Group& root, Polygon* p)
{
    using traits_type        = PolygonHDF5Traits;
    using h5_triangle_struct = traits_type::h5_triangle_struct;

    // ----------------------------------------------------------------------
    // load edge_lengths before assigning triangles
    Real3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);

    p->reset(edge_lengths);

    H5::DataSet triangle_dset(root.openDataSet("triangles"));
    const unsigned int num_triangles(
        triangle_dset.getSpace().getSimpleExtentNpoints());
    std::unique_ptr<h5_triangle_struct[]> h5_triangle_table(
        new h5_triangle_struct[num_triangles]);

    triangle_dset.read(h5_triangle_table.get(),
                       traits_type::get_triangle_comp_type());
    triangle_dset.close();

    std::vector<Triangle> triangles;
    triangles.reserve(num_triangles);
    for(std::size_t i=0; i<num_triangles; ++i)
    {
        Real3 p1, p2, p3;
        p1[0] = h5_triangle_table[i].p1x;
        p1[1] = h5_triangle_table[i].p1y;
        p1[2] = h5_triangle_table[i].p1z;
        p2[0] = h5_triangle_table[i].p2x;
        p2[1] = h5_triangle_table[i].p2y;
        p2[2] = h5_triangle_table[i].p2z;
        p3[0] = h5_triangle_table[i].p3x;
        p3[1] = h5_triangle_table[i].p3y;
        p3[2] = h5_triangle_table[i].p3z;
        triangles.push_back(Triangle(p1, p2, p3));
    }
    p->assign(triangles);
    return;
}

} // ecell4
