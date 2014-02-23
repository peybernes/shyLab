/// \file Serialization header. For now, only boost property tree is
/// used. We dump classes in structs in XML format.
#ifndef SERIALIZE_HPP
#define SERIALIZE_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

using boost::property_tree::ptree;

#endif
