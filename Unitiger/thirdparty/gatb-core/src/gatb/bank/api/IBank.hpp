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

/** \file IBank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK_IBANK_HPP_
#define _GATB_CORE_BANK_IBANK_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/bank/api/Sequence.hpp>

/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writting
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/********************************************************************************/

/** \brief Interface for what we need to read genomic databases.
 *
 * The IBank interface provides means to:
 *  - read sequences from some container.
 *  - insert sequences into some container.
 *
 * Typical implementation of this interface is a FASTA bank parser.
 *
 * The key point here is that clients should use this interface instead of specific
 * implementations, wherever they need to get nucleotides sequences from some container.
 *
 * By doing this, such clients will support many bank formats (according to the fact that
 * a IBank implementation provides such a service). They will also support "fake" banks, for
 * instance random generated sequences or any kind of synthetic data.
 */
class IBank : public tools::collections::Iterable<Sequence>, public tools::collections::Bag<Sequence>
{
public:

    /** Get an unique identifier for the bank (could be an URI for instance).
     * \return the identifier */
    virtual std::string getId () = 0;

    /** \copydoc tools::collections::Iterable::iterator */
    virtual tools::dp::Iterator<Sequence>* iterator () = 0;

    /** \copydoc tools::collections::Bag */
    virtual void insert (const Sequence& item) = 0;

    /** Return the size of the bank (comments + data)
     *
     * The returned value may be an approximation in some case. For instance, if we use
     * a zipped bank, an implementation may be not able to give accurate answer to the
     * size of the original file.
     *
     * \return the bank size.*/
    virtual u_int64_t getSize () = 0;

    /** Give an estimation of sequences information in the bank:
     *      - sequences number
     *      - sequences size (in bytes)
     *      - max size size (in bytes)
     * \return the sequences number estimation. */
    virtual void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize) = 0;

    /** Shortcut to 'estimate' method.
     * \return estimation of the number of sequences */
    virtual int64_t estimateNbItems () = 0;

    /** Shortcut to 'estimate' method.
     * \return estimation of the size of sequences */
    virtual u_int64_t estimateSequencesSize () = 0;

    /** \return the number of sequences read from the bank for computing estimated information */
    virtual u_int64_t getEstimateThreshold () = 0;

    /** Set the number of sequences read from the bank for computing estimated information
     * \param[in] nbSeq : the number of sequences to be read.*/
    virtual void setEstimateThreshold (u_int64_t nbSeq) = 0;
};

/********************************************************************************/

/** \brief Factory for IBank.
 *
 * This interface provides a factory method that builds a IBank* instance given some
 * identifier.
 *
 * Such an identifier can be an uri (FASTA banks for instance), or any mechanism allowing
 * to retrieve enough information for creating instances of a specific IBank implementation.
 */
class IBankFactory : public system::SmartPointer
{
public:

    /** Create an instance of IBank for a given uri.
     * \param[in] uri : the uri used for create the bank
     * \return the IBank instance. */
    virtual IBank* createBank (const std::string& uri) = 0;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IBANK_HPP_ */
