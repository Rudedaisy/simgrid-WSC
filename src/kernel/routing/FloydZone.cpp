/* Copyright (c) 2009-2022. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include <omp.h>

#include <climits>
#include <iostream>
#include <simgrid/kernel/routing/FloydZone.hpp>
#include <simgrid/kernel/routing/NetPoint.hpp>
#include <xbt/string.hpp>

#include "src/kernel/resource/NetworkModel.hpp"

#ifdef _COMPRESS
#include "bxzstr.hpp"
#endif

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(ker_routing_floyd, ker_routing, "Kernel Floyd Routing");

namespace simgrid {
namespace kernel::routing {

void FloydZone::init_tables(unsigned int table_size)
{
  if (link_table_.size() != table_size) {
    /* Resize and initialize Cost, Predecessor and Link tables */
    cost_table_.resize(table_size);
    link_table_.resize(table_size);
    predecessor_table_.resize(table_size);
    for (auto& cost : cost_table_)
      cost.resize(table_size, ULONG_MAX); /* link cost from host to host */
    for (auto& link : link_table_)
      link.resize(table_size); /* actual link between src and dst */
    for (auto& predecessor : predecessor_table_)
      predecessor.resize(table_size, -1); /* predecessor host numbers */
  }
}

void FloydZone::get_local_route(const NetPoint* src, const NetPoint* dst, Route* route, double* lat)
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

void FloydZone::add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                          const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical)
{
  /* set the size of table routing */
  unsigned int table_size = get_table_size();
  init_tables(table_size);

  add_route_check_params(src, dst, gw_src, gw_dst, link_list, symmetrical);

  /* Check that the route does not already exist */
  if (gw_dst && gw_src) // netzone route (to adapt the error message, if any)
    xbt_assert(nullptr == link_table_[src->id()][dst->id()],
               "The route between %s@%s and %s@%s already exists (Rq: "
               "routes are symmetrical by default).",
               src->get_cname(), gw_src->get_cname(), dst->get_cname(), gw_dst->get_cname());
  else
    xbt_assert(nullptr == link_table_[src->id()][dst->id()],
               "The route between %s and %s already exists (Rq: routes are "
               "symmetrical by default).",
               src->get_cname(), dst->get_cname());

  link_table_[src->id()][dst->id()] = std::unique_ptr<Route>(
      new_extended_route(get_hierarchy(), gw_src, gw_dst, get_link_list_impl(link_list, false), true));
  predecessor_table_[src->id()][dst->id()] = src->id();
  cost_table_[src->id()][dst->id()]        = link_table_[src->id()][dst->id()]->link_list_.size();

  if (symmetrical) {
    if (gw_dst && gw_src) // netzone route (to adapt the error message, if any)
      xbt_assert(nullptr == link_table_[dst->id()][src->id()],
                 "The route between %s@%s and %s@%s already exists. You "
                 "should not declare the reverse path as symmetrical.",
                 dst->get_cname(), gw_dst->get_cname(), src->get_cname(), gw_src->get_cname());
    else
      xbt_assert(nullptr == link_table_[dst->id()][src->id()],
                 "The route between %s and %s already exists. You should "
                 "not declare the reverse path as symmetrical.",
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
    cost_table_[dst->id()][src->id()] =
        link_table_[dst->id()][src->id()]->link_list_.size(); /* count of links, old model assume 1 */
  }
}

// aux functions for floyd routing
unsigned long* floyd_warshall_blocked_init(const int n, const int block_size, int* n_blocked)
{
  int n_oversized;
  int block_remainder = n % block_size;
  if (block_remainder == 0) {
    n_oversized = n;
  } else {
    n_oversized = n + block_size - block_remainder;
  }
  *n_blocked         = n_oversized;
  unsigned long* out = new unsigned long[n_oversized * n_oversized];
  return out;
}

inline void floyd_warshall_in_place(unsigned long* C, const unsigned long* A, const unsigned long* B,
                                    signed long* aux_C, signed long* aux_B, const int ori_n, const int b, const int n,
                                    const int k_base, const int i_base, const int j_base, const int host_id,
                                    const int io_router_id)
{
  for (int k = 0; k < b; k++) {
    if (k + k_base >= ori_n)
      break;
    int kth = k * n;
    for (int i = 0; i < b; i++) {
      if (i + i_base >= ori_n)
        break;
      for (int j = 0; j < b; j++) {
        if (j + j_base >= ori_n)
          break;
        if (k + k_base == io_router_id && i + i_base != host_id && j + j_base != host_id &&
            i + i_base != io_router_id && j + j_base != io_router_id)
          continue;
        unsigned long sum = A[i * n + k] + B[kth + j];
        if (A[i * n + k] < ULONG_MAX && B[kth + j] < ULONG_MAX && (C[i * n + j] > sum || C[i * n + j] == ULONG_MAX)) {
          C[i * n + j]     = sum;
          aux_C[i * n + j] = aux_B[kth + j];
        }
      }
    }
  }
}

void floyd_warshall_blocked(unsigned long* output, signed long* aux, const int ori_n, const int n, const int b,
                            const int io_router_id, const int host_id)
{
  // for now, assume b divides n
  const int blocks = n / b;

  // note that [i][j] == [i * input_width * block_width + j * block_width]
  for (int k = 0; k < blocks; k++) {
    floyd_warshall_in_place(&output[k * b * n + k * b], &output[k * b * n + k * b], &output[k * b * n + k * b],
                            &aux[k * b * n + k * b], &aux[k * b * n + k * b], ori_n, b, n, k * b, k * b, k * b, host_id,
                            io_router_id);
#pragma omp parallel for
    for (int j = 0; j < blocks; j++) {
      if (j == k)
        continue;
      floyd_warshall_in_place(&output[k * b * n + j * b], &output[k * b * n + k * b], &output[k * b * n + j * b],
                              &aux[k * b * n + j * b], &aux[k * b * n + j * b], ori_n, b, n, k * b, k * b, j * b,
                              host_id, io_router_id);
    }
#pragma omp parallel for
    for (int i = 0; i < blocks; i++) {
      if (i == k)
        continue;
      floyd_warshall_in_place(&output[i * b * n + k * b], &output[i * b * n + k * b], &output[k * b * n + k * b],
                              &aux[i * b * n + k * b], &aux[k * b * n + k * b], ori_n, b, n, k * b, i * b, k * b,
                              host_id, io_router_id);
      for (int j = 0; j < blocks; j++) {
        if (j == k)
          continue;
        floyd_warshall_in_place(&output[i * b * n + j * b], &output[i * b * n + k * b], &output[k * b * n + j * b],
                                &aux[i * b * n + j * b], &aux[k * b * n + j * b], ori_n, b, n, k * b, i * b, j * b,
                                host_id, io_router_id);
      }
    }
  }
}

void FloydZone::do_seal()
{
  /* Load checkpoint if exists, then skip executing this function */
  auto ckpt_path = this->get_property("ckpt_path");
#ifdef _COMPRESS
  bxz::ifstream ifs(ckpt_path);
#else
  /* Load checkpoint if exists, then skip executing this function */
  std::ifstream ifs(ckpt_path, std::ios::binary);
#endif
  if (ifs.good()) {
    boost::archive::text_iarchive ia(ifs);
    ia >> *this;
    return;
  }

  std::cout << "Routing Checkpoint not found! Generating now." << std::endl;

  /* set the size of table routing */
  unsigned int table_size = get_table_size();
  init_tables(table_size);
  unsigned host_id      = 0;
  unsigned io_router_id = table_size - 1;

  /* Add the loopback if needed */
  if (get_network_model()->loopback_ && get_hierarchy() == RoutingMode::base) {
    for (unsigned int i = 0; i < table_size; i++) {
      auto& route = link_table_[i][i];
      if (not route) {
        route.reset(new Route());
        route->link_list_.push_back(get_network_model()->loopback_.get());
        predecessor_table_[i][i] = i;
        cost_table_[i][i]        = 1;
      }
    }
  }

  // Blocked Parallel Floyd
  // Reference: https://moorejs.github.io/APSP-in-parallel/
  int blocked_table_size = 0;
  int block_size         = 256;
  unsigned long* _ct     = floyd_warshall_blocked_init(table_size, block_size, &blocked_table_size);
  signed long* _pret     = (signed long*)floyd_warshall_blocked_init(table_size, block_size, &blocked_table_size);

  for (int i = 0; i < table_size; i++) {
    memcpy(_ct + i * blocked_table_size, cost_table_[i].data(), sizeof(unsigned long) * table_size);
    memcpy(_pret + i * blocked_table_size, predecessor_table_[i].data(), sizeof(unsigned long) * table_size);
  }

  floyd_warshall_blocked(_ct, _pret, table_size, blocked_table_size, block_size, io_router_id, host_id);

  for (int i = 0; i < table_size; i++) {
    memcpy(cost_table_[i].data(), _ct + i * blocked_table_size, sizeof(unsigned long) * table_size);
    memcpy(predecessor_table_[i].data(), _pret + i * blocked_table_size, sizeof(unsigned long) * table_size);
  }

  free(_ct);
  free(_pret);
  // TODO: Modify the checkpoint file's name
#ifdef _COMPRESS
  bxz::ofstream ofs(ckpt_path, bxz::lzma, 5);
  if (!ofs.good()) {
    std::cout << "Invalid checkpoint path: " << ckpt_path << ". Use default." << std::endl;
    bxz::ofstream ofs("FloydRoutingCheckpoint.txt", bxz::lzma, 5);
  }
#else
  std::ofstream ofs("FloydRoutingCheckpoint.txt");
  if (!ofs.good()) {
    std::cout << "Invalid checkpoint path: " << ckpt_path << ". Use default." << std::endl;
    std::ofstream ofs("FloydRoutingCheckpoint.txt", bxz::lzma, 5);
  }
#endif
  boost::archive::text_oarchive oa(ofs);
  oa << *this;
}
} // namespace kernel::routing

namespace s4u {
NetZone* create_floyd_zone(const std::string& name)
{
  return (new kernel::routing::FloydZone(name))->get_iface();
}
} // namespace s4u

} // namespace simgrid
