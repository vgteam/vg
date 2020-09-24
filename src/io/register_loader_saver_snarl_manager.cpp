/**
 * \file register_loader_saver_snarl_manager.cpp
 * Defines IO for a SnarlManager from stream files.
 */

#include <vg/io/registry.hpp>
#include <vg/io/protobuf_iterator.hpp>
#include "register_loader_saver_snarl_manager.hpp"

#include "../snarls.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_snarl_manager() {
    Registry::register_loader_saver<SnarlManager>(vector<string>{"SNARL", ""}, [](const message_sender_function_t& for_each_message) -> void* {
       SnarlManager* manager = new SnarlManager([&](const function<void(Snarl&)>& consume_snarl) {
            // Call the source function with a function that deserializes each message and feeds it to the SnarlManager.
            for_each_message([&](const string& serialized_snarl) {
                // Parse the message to a Snarl
                Snarl s;
                if (!ProtobufIterator<Snarl>::parse_from_string(s, serialized_snarl)) {
                    // Handle bad graphs.
                    // TODO: make this an exception if we ever want to be allowed to continue from this.
                    cerr << "error[register_loader_saver_snarl_manager]: invalid Snarl message" << endl;
                    exit(1);
                }
                
                // Feed the Snarl to the SnarlManager.
                // TODO: Is all this callback-and-callforth better than just a loop somewhere?
                // TODO: Unify with the VG load logic to avoid duplicating deserialization
                consume_snarl(s);
            });
        });
        
        return (void*) manager;
    }, [](const void* manager_void, const message_consumer_function_t& send_message) {
        // Cast to SnarlManager and serialize to the consumer.
        assert(manager_void != nullptr);
        
        const SnarlManager* manager = (const SnarlManager*) manager_void;
        
        manager->for_each_snarl_preorder([&](const Snarl* snarl) {
            // Serialize and emit each snarl
            string s;
            snarl->SerializeToString(&s);
            send_message(s);
        });
    });
}

}

}

