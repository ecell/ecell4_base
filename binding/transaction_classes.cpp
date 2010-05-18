#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Transaction.hpp"
#include "binding_common.hpp"

namespace binding {

void register_transaction_classes()
{
    register_transaction_class<Transaction, ParticleContainer>("Transaction");
    register_transaction_impl_class<TransactionImpl, Transaction>("TransactionImpl");
}

} // namesapce binding
