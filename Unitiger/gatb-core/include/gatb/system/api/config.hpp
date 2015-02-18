/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#define INT128_FOUND            1

#define STR_LIBRARY_VERSION     "1.0.5"
#define STR_COMPILATION_DATE    "2014-09-09 16:01:54"
#define STR_COMPILATION_FLAGS   " -O3  -DINT128_FOUND  -DWITH_LAMBDA_EXPRESSION  -std=c++0x -DBOOST_NO_CXX11_RVALUE_REFERENCES=1 -Wno-invalid-offsetof"
#define STR_COMPILER            "/usr/bin/cc  (4.2.0)"

#define STR_OPERATING_SYSTEM    "Darwin-11.4.2"

#define KSIZE_1  32
#define KSIZE_2  64
#define KSIZE_3  96   
#define KSIZE_4  128

#define PREC_1  ((KSIZE_1+31)/32)
#define PREC_2  ((KSIZE_2+31)/32)
#define PREC_3  ((KSIZE_3+31)/32)
#define PREC_4  ((KSIZE_4+31)/32)
