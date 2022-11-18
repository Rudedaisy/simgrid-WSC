/* Copyright (c) 2009-2022. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include <simgrid/kernel/routing/WSCZone.hpp>
#include <simgrid/kernel/routing/NetPoint.hpp>
#include <xbt/string.hpp>
//#include <iostream>
#include <cstdlib> // for std::abs()
//#include <boost/stacktrace.hpp>
//#include "src/kernel/resource/StandardLinkImpl.hpp"
#include "src/kernel/resource/NetworkModel.hpp"
#include <omp.h> // speed up the seal() function

#include <climits>

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(ker_routing_wsc, ker_routing, "Kernel WSC Routing");

namespace simgrid {
namespace kernel {
namespace routing {

void WSCZone::init_tables(unsigned int table_size)
{
  // Table size = #compute_nodes + 1
  if (link_table_.size() != table_size) {
    /* Resize and initialize Cost, Predecessor and Link tables */
    link_table_.resize(table_size);
    predecessor_table_.resize(table_size);
    for (auto& link : link_table_)
      link.resize(table_size); /* actual link between src and dst */
    for (auto& predecessor : predecessor_table_)
      predecessor.resize(table_size, -1); /* predecessor host numbers */
  }
}

void WSCZone::get_local_route(const NetPoint* src, const NetPoint* dst, Route* route, double* lat)
{
  get_route_check_params(src, dst);

  /* create a result route */
  std::vector<Route*> route_stack;
  unsigned long cur = dst->id();
  do {
    long pred = predecessor_table_[src->id()][cur];
    if (pred == -1)
      throw std::invalid_argument(xbt::string_printf("No route from '%s' to '%s'", src->get_cname(), dst->get_cname()));
    route_stack.push_back(link_table_[pred][cur].get());
    cur = pred;
  } while (cur != src->id());

  if (get_hierarchy() == RoutingMode::recursive) {
    route->gw_src_ = route_stack.back()->gw_src_;
    route->gw_dst_ = route_stack.front()->gw_dst_;
  }

  const NetPoint* prev_dst_gw = nullptr;
  while (not route_stack.empty()) {
    const Route* e_route = route_stack.back();
    route_stack.pop_back();
    if (get_hierarchy() == RoutingMode::recursive && prev_dst_gw != nullptr &&
        prev_dst_gw->get_cname() != e_route->gw_src_->get_cname()) {
      get_global_route(prev_dst_gw, e_route->gw_src_, route->link_list_, lat);
    }
    
    add_link_latency(route->link_list_, e_route->link_list_, lat);

    prev_dst_gw = e_route->gw_dst_;
  }
}

void WSCZone::add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                          const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical)
{
  
  /* set the size of table routing */
  unsigned int table_size = get_table_size();
  init_tables(table_size);

  add_route_check_params(src, dst, gw_src, gw_dst, link_list, symmetrical);

  /* Check that the route does not already exist */
  if (gw_dst && gw_src) // netzone route (to adapt the error message, if any)
    xbt_assert(nullptr == link_table_[src->id()][dst->id()],
               "The route between %s@%s and %s@%s already exists (Rq: routes are symmetrical by default).",
               src->get_cname(), gw_src->get_cname(), dst->get_cname(), gw_dst->get_cname());
  else
    xbt_assert(nullptr == link_table_[src->id()][dst->id()],
               "The route between %s and %s already exists (Rq: routes are symmetrical by default).", src->get_cname(),
               dst->get_cname());

  link_table_[src->id()][dst->id()] = std::unique_ptr<Route>(
      new_extended_route(get_hierarchy(), gw_src, gw_dst, get_link_list_impl(link_list, false), true));
  predecessor_table_[src->id()][dst->id()] = src->id();

  if (symmetrical) {
    if (gw_dst && gw_src) // netzone route (to adapt the error message, if any)
      xbt_assert(
          nullptr == link_table_[dst->id()][src->id()],
          "The route between %s@%s and %s@%s already exists. You should not declare the reverse path as symmetrical.",
          dst->get_cname(), gw_dst->get_cname(), src->get_cname(), gw_src->get_cname());
    else
      xbt_assert(nullptr == link_table_[dst->id()][src->id()],
                 "The route between %s and %s already exists. You should not declare the reverse path as symmetrical.",
                 dst->get_cname(), src->get_cname());

    if (gw_dst && gw_src) {
      NetPoint* gw_tmp = gw_src;
      gw_src           = gw_dst;
      gw_dst           = gw_tmp;
    }

    if (not gw_src || not gw_dst)
      XBT_DEBUG("Load Route from \"%s\" to \"%s\"", dst->get_cname(), src->get_cname());
    else
      XBT_DEBUG("Load NetzoneRoute from \"%s(%s)\" to \"%s(%s)\"", dst->get_cname(), gw_src->get_cname(),
                src->get_cname(), gw_dst->get_cname());

    link_table_[dst->id()][src->id()] = std::unique_ptr<Route>(
        new_extended_route(get_hierarchy(), gw_src, gw_dst, get_link_list_impl(link_list, true), false));
    predecessor_table_[dst->id()][src->id()] = dst->id();
  }

  /* Roundabout way to infer WSC dimensions. Find X dimension by the order of adding routes. CAUTION: unstable to loop order of gen_config */
  if(!dimensions_inferred && src->id() == contiguous_x && dst->id() == contiguous_x+1) contiguous_x++;
  
  //std::cout << "Linking " << src->get_cname() << "(" << src->id() << ") to " << dst->get_cname() << "(" << dst->id() << ")" << std::endl;
}

void WSCZone::do_seal()
{
  
  /* set the size of table routing */
  unsigned int table_size = get_table_size();
  init_tables(table_size);

  if(!dimensions_inferred) {
    wsc_x = contiguous_x;
    wsc_y = (table_size-2) / wsc_x;
    //std::cout << "Inferred wafer dims " << wsc_x << "x" << wsc_y << std::endl;
    dimensions_inferred = true;
  }

  /* Add the loopback if needed */
  if (get_network_model()->loopback_ && get_hierarchy() == RoutingMode::base) {
    for (unsigned int i = 0; i < table_size; i++) {
      auto& route = link_table_[i][i];
      if (not route) {
        route.reset(new Route());
        route->link_list_.push_back(get_network_model()->loopback_.get());
        predecessor_table_[i][i] = i;
      }
    }
  }
  
  long int host1_id = 0;
  long int io_router_id = table_size-1;
  long int perim_dists[4] = {0}; // used to find which perimeter is closest to dst_pe. Format: [N,S,W,E]
  long int src_pe;
  long int dst_pe;
  /* Connect all nodes via predecessor_table_ */
  /* Special nodes: 0=host1, table_size-1=io_router */
  int direction_x = 1;
  int direction_y = 1;
  predecessor_table_[0][io_router_id] = 0;
#pragma omp parallel for collapse(2) num_threads(OPENMP_THREADS) shared(predecessor_table_) private(src_pe,dst_pe,direction_x,direction_y,perim_dists)
  for (long int src_y = 0; src_y < wsc_y; src_y++) {
    for (long int src_x = 0; src_x < wsc_x; src_x++) {
      src_pe = src_y*wsc_x + src_x + 1; // +1 for host1 offset
      /* First: internal nodes */
      for (long int dst_y = src_y; dst_y < wsc_y; dst_y++) {
	direction_y = dst_y==src_y ? 0 : (dst_y-src_y) / std::abs(dst_y-src_y);
	for (long int dst_x = (dst_y==src_y ? src_x+1 : 0); dst_x < wsc_x; dst_x++) {
	  dst_pe = dst_y*wsc_x + dst_x + 1; // +1 for host1 offset
	  direction_x = dst_x==src_x ? 0 : (dst_x-src_x) / std::abs(dst_x-src_x);

	  if (direction_x != 0) {
	    predecessor_table_[src_pe][dst_pe] = dst_pe - direction_x;
	    predecessor_table_[dst_pe][src_pe] = src_pe + direction_x;
	  } else {
	    predecessor_table_[src_pe][dst_pe] = dst_pe - (direction_y*wsc_x);
	    predecessor_table_[dst_pe][src_pe] = src_pe + (direction_y*wsc_x);
	  }
	}
      }

      /* Next: host-to-PEs */
      dst_pe = src_pe;
      long int dst_y = src_y;
      long int dst_x = src_x;
      perim_dists[0] = dst_y;
      perim_dists[1] = wsc_y-dst_y-1;
      perim_dists[2] = dst_x;
      perim_dists[3] = wsc_x-dst_x-1;

      // Different routings depending on distance from the perimeters
      if (perim_dists[0]<=perim_dists[1] && perim_dists[0]<=perim_dists[2] && perim_dists[0]<=perim_dists[3]) {
        // N
        predecessor_table_[host1_id][dst_pe] = dst_y==0 ? io_router_id : dst_pe - wsc_x;
        predecessor_table_[dst_pe][host1_id] = io_router_id;
        predecessor_table_[io_router_id][dst_pe] = dst_y==0 ? io_router_id : dst_pe - wsc_x;
        predecessor_table_[dst_pe][io_router_id] = dst_y==0 ? dst_pe : dst_x+1; //dst_pe % wsc_x; // row-0, same column                                                                       
      } else if (perim_dists[1]<=perim_dists[0] && perim_dists[1]<=perim_dists[2] && perim_dists[1]<=perim_dists[3]) {
        // S
        predecessor_table_[host1_id][dst_pe] = dst_y==wsc_y-1 ? io_router_id : dst_pe + wsc_x;
        predecessor_table_[dst_pe][host1_id] = io_router_id;
        predecessor_table_[io_router_id][dst_pe] = dst_y==wsc_y-1 ? io_router_id : dst_pe + wsc_x;
        predecessor_table_[dst_pe][io_router_id] = dst_y==wsc_y-1 ? dst_pe : table_size-2 - wsc_x + dst_x+1;
      } else if (perim_dists[2]<=perim_dists[0] && perim_dists[2]<=perim_dists[1] && perim_dists[2]<=perim_dists[3]) {
        // W
        predecessor_table_[host1_id][dst_pe] = dst_x==0 ? io_router_id : dst_pe - 1;
        predecessor_table_[dst_pe][host1_id] = io_router_id;
        predecessor_table_[io_router_id][dst_pe] = dst_x==0 ? io_router_id : dst_pe - 1;
        predecessor_table_[dst_pe][io_router_id] = dst_x==0 ? dst_pe : dst_y*wsc_x+1;
      } else {
        // E
        predecessor_table_[host1_id][dst_pe] = dst_x==wsc_x-1 ? io_router_id : dst_pe + 1;
        predecessor_table_[dst_pe][host1_id] = io_router_id;
        predecessor_table_[io_router_id][dst_pe] = dst_x==wsc_x-1 ? io_router_id : dst_pe + 1;
        predecessor_table_[dst_pe][io_router_id] = dst_x==wsc_x-1 ? dst_pe : dst_y*wsc_x + wsc_x;
      }
    }
  }

} // do_seal()
} // namespace routing
} // namespace kernel

namespace s4u {
NetZone* create_wsc_zone(const std::string& name)
{
  return (new kernel::routing::WSCZone(name))->get_iface();
}
} // namespace s4u

} // namespace simgrid
