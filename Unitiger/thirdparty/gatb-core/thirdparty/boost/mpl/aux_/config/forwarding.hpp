
#ifndef BOOST_MPL_AUX_CONFIG_FORWARDING_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_FORWARDING_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: forwarding.hpp,v 1.1 2014/12/19 07:36:31 lmsalmel Exp $
// $Date: 2014/12/19 07:36:31 $
// $Revision: 1.1 $

#include <boost/mpl/aux_/config/workaround.hpp>

#if    !defined(BOOST_MPL_CFG_NO_NESTED_FORWARDING) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE) \
    && BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x610))

#   define BOOST_MPL_CFG_NO_NESTED_FORWARDING

#endif

#endif // BOOST_MPL_AUX_CONFIG_FORWARDING_HPP_INCLUDED
