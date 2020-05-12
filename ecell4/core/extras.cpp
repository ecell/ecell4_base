#include "extras.hpp"
#include "exceptions.hpp"
#include <iostream>
#include <sstream>
#include <regex>


namespace ecell4
{

namespace extras
{

Shape::dimension_kind
get_dimension_from_model(const Species& species, const std::shared_ptr<Model>& model)
{
    const Shape::dimension_kind DEFAULT_DIMENSION(Shape::THREE);

    if (species.serial().empty())
        return DEFAULT_DIMENSION;

    if (!model->has_species_attribute(species))
    {
        // std::stringstream ss;
        // ss << "The model has no attribute for Specis(\"" << species.serial() << "\")";
        // throw NotFound(ss.str());
        return DEFAULT_DIMENSION;
    }

    const Species& attribute(model->apply_species_attributes(species));

    if (attribute.has_attribute("dimension"))
    {
        switch (attribute.get_attribute_as<Integer>("dimension"))
        {
            case 1: return Shape::ONE;
            case 2: return Shape::TWO;
            case 3: return Shape::THREE;
        }
    }

    if (attribute.has_attribute("location"))
    {
        return get_dimension_from_model(
            Species(attribute.get_attribute_as<std::string>("location")),
            model
        );
    }

    return DEFAULT_DIMENSION;
}

#ifdef WITH_HDF5
void save_version_information(H5::H5Location* root, const std::string& version)
{
    if (version.size() > 32)
    {
        throw IllegalArgument("Version info must be shorter than 32 characters.");
    }
    using namespace H5;
    std::unique_ptr<DataSet> dataset(
        new DataSet(root->createDataSet("version", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))));
    dataset->write(version.c_str(), dataset->getDataType());
}

std::string load_version_information(const H5::H5Location& root)
{
    using namespace H5;
    char buf[32];
    const DataSet dataset(DataSet(root.openDataSet("version")));
    dataset.read(buf, dataset.getDataType());
    return std::string(buf);
}
#endif

std::string load_version_information(const std::string& filename)
{
#ifdef WITH_HDF5
    using namespace H5;
    std::unique_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    return load_version_information(*fin);
#else
    return "";
#endif
}

template <typename T>
T mystoi(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    T retval;
    ss >> retval;
    return retval;
}

std::pair<VersionInformation::prerelease_type, unsigned int> parse_prerelease(const std::string& prestr)
{
    if (prestr.size() == 0)
    {
        return std::make_pair(VersionInformation::FINAL, 0);
    }
    else if (prestr[0] == 'a')
    {
        return std::make_pair(VersionInformation::ALPHA, mystoi<unsigned int>(prestr.substr(1)));
    }
    else if (prestr[0] == 'b')
    {
        return std::make_pair(VersionInformation::BETA, mystoi<unsigned int>(prestr.substr(1)));
    }
    else if (prestr[0] == 'c')
    {
        return std::make_pair(VersionInformation::RC, mystoi<unsigned int>(prestr.substr(1)));
    }
    else if (prestr.size() >= 2 && prestr[0] == 'r' && prestr[1] == 'c')
    {
        return std::make_pair(VersionInformation::RC, mystoi<unsigned int>(prestr.substr(2)));
    }
    else
    {
        throw NotSupported("Unknown pre-release was given.");
        return std::make_pair(VersionInformation::NONE, 0);
    }
}

VersionInformation parse_version_information(const std::string& version)
{
    using namespace std;

    regex reg("^([^-\\.]+-[^-\\.]+-)([0123456789]+)\\.([0123456789]+)(\\.[0123456789]+|)(a[0123456789]+|b[0123456789]+|rc[0123456789]+|c[0123456789]+|)(\\.dev[0123456789]+|)$");
    smatch result;
    if (!regex_match(version, result, reg))
    {
        throw std::invalid_argument(
            "a wrong version information was given [" + version + "]"); //XXX:
    }

    const std::string header = result.str(1);
    const unsigned int majorno = mystoi<unsigned int>(result.str(2));
    const unsigned int minorno = mystoi<unsigned int>(result.str(3));
    const unsigned int patchno = (result.str(4).size() > 1 ? mystoi<unsigned int>(result.str(4).substr(1)) : 0);
    const std::pair<VersionInformation::prerelease_type, unsigned int> pre = parse_prerelease(result.str(5));
    const int devno = (result.str(6).size() > 4 ? mystoi<int>(result.str(6).substr(4)) : -1);

    return VersionInformation(header, majorno, minorno, patchno, pre.first, pre.second, devno);

}

bool check_version_information(const std::string& version, const std::string& required)
{
    const VersionInformation vinfo1(parse_version_information(version));
    const VersionInformation vinfo2(parse_version_information(required));

    if (vinfo1.header != vinfo2.header)
    {
        return false;
    }
    else if (vinfo1.majorno != vinfo2.majorno)
    {
        return (vinfo1.majorno > vinfo2.majorno);
    }
    else if (vinfo1.minorno != vinfo2.minorno)
    {
        return (vinfo1.minorno > vinfo2.minorno);
    }
    else if (vinfo1.patchno != vinfo2.patchno)
    {
        return (vinfo1.patchno > vinfo2.patchno);
    }
    else if (vinfo1.pre != vinfo2.pre)
    {
        return (vinfo1.pre > vinfo2.pre);
    }
    else if (vinfo1.preno != vinfo2.preno)
    {
        return (vinfo1.preno > vinfo2.preno);
    }
    return (vinfo1.devno == -1 || (vinfo2.devno != -1 && vinfo1.devno >= vinfo2.devno));
}

std::vector<std::vector<Real> > get_stoichiometry(
    const std::vector<Species>& species_list, const std::vector<ReactionRule>& reaction_rules)
{
    typedef std::unordered_map<Species, unsigned> species_map_type;

    species_map_type index_map;
    {
        unsigned i(0);
        for(std::vector<Species>::const_iterator it(species_list.begin());
            it != species_list.end(); it++, i++)
        {
            index_map[*it] = i;
        }
    }

    std::vector<std::vector<Real> > ret;
    {
        ret.resize(species_list.size());
        for (std::vector<std::vector<Real> >::iterator it(ret.begin());
             it != ret.end(); ++it)
        {
            (*it).resize(reaction_rules.size());
        }

        unsigned i(0);
        for (std::vector<ReactionRule>::const_iterator it(reaction_rules.begin());
             it != reaction_rules.end(); ++it, i++)
        {
            const ReactionRule& rr(*it);

            if (!rr.has_descriptor())
            {
                for (auto const& sp : rr.reactants())
                {
                    ret[index_map[sp]][i] -= 1.0;
                }
                for (auto const& sp : rr.products())
                {
                    ret[index_map[sp]][i] += 1.0;
                }
            }
            else
            {
                const std::shared_ptr<ReactionRuleDescriptor>& desc(rr.get_descriptor());
                {
                    if (rr.reactants().size() != desc->reactant_coefficients().size())
                    {
                        std::stringstream msg;
                        msg << "The number of reactant coefficients mismatches ("
                            << desc->reactant_coefficients().size()
                            << " != "
                            << rr.reactants().size()
                            << ").";
                        throw std::runtime_error(msg.str());
                    }

                    ReactionRule::reactant_container_type::const_iterator it1(rr.reactants().begin());
                    ReactionRuleDescriptor::coefficient_container_type::const_iterator it2(desc->reactant_coefficients().begin());
                    for (; it1 != rr.reactants().end(); ++it1, ++it2)
                    {
                        ret[index_map[(*it1)]][i] -= (*it2);
                    }
                }
                {
                    if (rr.products().size() != desc->product_coefficients().size())
                    {
                        std::stringstream msg;
                        msg << "The number of product coefficients mismatches ("
                            << desc->product_coefficients().size()
                            << " != "
                            << rr.products().size()
                            << ").";
                        throw std::runtime_error(msg.str());
                    }

                    ReactionRule::product_container_type::const_iterator it1(rr.products().begin());
                    ReactionRuleDescriptor::coefficient_container_type::const_iterator it2(desc->product_coefficients().begin());
                    for (; it1 != rr.products().end(); ++it1, ++it2)
                    {
                        ret[index_map[(*it1)]][i] += (*it2);
                    }
                }
            }
        }
    }
    return ret;
}

} // extras

} // ecell4
