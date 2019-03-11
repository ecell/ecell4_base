#include "extras.hpp"
#include "exceptions.hpp"
#include <iostream>
#include <sstream>

#if defined(HAVE_BOOST_REGEX)
#include <boost/regex.hpp>
#elif defined(WIN32_MSC)
#include <regex>
#else
#include <regex.h>
#endif /* HAVE_BOOST_REGEX */


namespace ecell4
{

namespace extras
{

Shape::dimension_kind
get_dimension_from_model(const Species& species, const boost::shared_ptr<Model>& model)
{
    const Shape::dimension_kind DEFAULT_DIMENSION(Shape::THREE);

    if (species.serial().empty())
        return DEFAULT_DIMENSION;

    if (!model->has_species_attribute(species))
    {
        std::stringstream ss;
        ss << "The model has no attribute for Specis(\"" << species.serial() << "\")";
        throw NotFound(ss.str());
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
    boost::scoped_ptr<DataSet> dataset(
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
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    return load_version_information(*fin);
#else
    return "";
#endif
}

template <typename T = int>
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
#if defined(HAVE_BOOST_REGEX) || defined(WIN32_MSC)
#if defined(HAVE_BOOST_REGEX)
    using namespace boost;
#else /* WIN32_MSC */
    using namespace std::tr1;
#endif /* HAVE_BOOST_REGEX */
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
#else /* regex.h */
    regex_t reg;
    int errcode = regcomp(
        &reg, "^([^-\\.]+-[^-\\.]+-)([0123456789]+)\\.([0123456789]+)(\\.[0123456789]+|)(a[0123456789]+|b[0123456789]+|rc[0123456789]+|c[0123456789]+|)(\\.dev[0123456789]+|)$",
        REG_EXTENDED);
    if (errcode != 0)
    {
        char errbuf[100];
        regerror(errcode, &reg, errbuf, sizeof(errbuf));
        regfree(&reg);
        std::cout << "regcompile error: " << errbuf << std::endl;
        throw IllegalState("regcompile error.");
    }

    regmatch_t match[7];
    errcode = regexec(&reg, version.c_str(), 7, match, 0);
    if (errcode != 0)
    {
        char errbuf[100];
        regerror(errcode, &reg, errbuf, sizeof(errbuf));
        regfree(&reg);
        std::cout << "regexec error: " << errbuf << std::endl;
        throw IllegalState("regexec error.");
    }

    const std::string header = version.substr(match[1].rm_so, match[1].rm_eo - match[1].rm_so);
    const unsigned int majorno = mystoi<unsigned int>(version.substr(match[2].rm_so, match[2].rm_eo - match[2].rm_so));
    const unsigned int minorno = mystoi<unsigned int>(version.substr(match[3].rm_so, match[3].rm_eo - match[3].rm_so));
    const unsigned int patchno = (match[4].rm_eo - match[4].rm_so > 0 ?
            mystoi<unsigned int>(version.substr(match[4].rm_so + 1, match[4].rm_eo - (match[4].rm_so + 1))) : 0);
    const std::pair<VersionInformation::prerelease_type, unsigned int> pre = parse_prerelease(
            version.substr(match[5].rm_so, match[5].rm_eo - match[5].rm_so));
    const int devno = (match[6].rm_eo - match[6].rm_so > 0 ?
            mystoi<int>(version.substr(match[6].rm_so + 4, match[6].rm_eo - (match[6].rm_so + 4))) : -1);

    regfree(&reg);

    return VersionInformation(header, majorno, minorno, patchno, pre.first, pre.second, devno);
#endif /* HAVE_BOOST_REGEX */
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

} // extras

} // ecell4
