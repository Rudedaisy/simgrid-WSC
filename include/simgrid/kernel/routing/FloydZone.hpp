/* Copyright (c) 2013-2022. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SURF_ROUTING_FLOYD_HPP_
#define SURF_ROUTING_FLOYD_HPP_

#include <simgrid/kernel/routing/RoutedZone.hpp>
#include <fstream>
// include headers that implement a archive in simple text format
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace simgrid {
namespace kernel {
namespace routing {

/** @ingroup ROUTING_API
 *  @brief NetZone with an explicit routing computed at initialization with Floyd-Warshal
 *
 *  The path between components is computed at creation time from every one-hop links,
 *  using the Floyd-Warshal algorithm.
 *
 *  This result in rather small platform file, slow initialization time,  and intermediate memory requirements
 *  (somewhere between the one of @{DijkstraZone} and the one of @{FullZone}).
 */
class XBT_PRIVATE FloydZone : public RoutedZone {
  /* members for serialization of tables */
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
    ar & predecessor_table_;
    ar & cost_table_;
    //ar & link_table_;
  }
  
  /* vars to compute the Floyd algorithm. */
  std::vector<std::vector<long>> predecessor_table_;
  std::vector<std::vector<unsigned long>> cost_table_;
  std::vector<std::vector<std::unique_ptr<Route>>> link_table_;

  void init_tables(unsigned int table_size);
  void do_seal() override;

public:
  using RoutedZone::RoutedZone;
  FloydZone(const FloydZone&) = delete;
  FloydZone& operator=(const FloydZone&) = delete;

  void get_local_route(const NetPoint* src, const NetPoint* dst, Route* into, double* latency) override;
  void add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                 const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical) override;
};
} // namespace routing
} // namespace kernel
} // namespace simgrid

#endif /* SURF_ROUTING_FLOYD_HPP_ */
