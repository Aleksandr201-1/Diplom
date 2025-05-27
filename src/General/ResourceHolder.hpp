#ifndef RESOURCE_HOLDER_HPP
#define RESOURCE_HOLDER_HPP

#include <string>
#include <memory>

template <class Resource, class Identifier>
class ResourceHolder {
    public:
        ResourceHolder ();
        ~ResourceHolder ();
        void load (Identifier id, std::unique_ptr<Resource> &resource);
        void load (Identifier id, const std::string& filename);
        template <class Parameter>
        void load (Identifier id, const std::string& filename, const Parameter& secondParam);
        Resource &get (Identifier id);
        const Resource &get (Identifier id) const;

        Resource &operator[] (Identifier id);
        const Resource &operator[] (Identifier id) const;
    private:
        std::map<Identifier, std::unique_ptr<Resource>> resourceMap;
};

template <class Resource, class Identifier>
ResourceHolder<Resource, Identifier>::ResourceHolder () {}

template <class Resource, class Identifier>
ResourceHolder<Resource, Identifier>::~ResourceHolder () {}

template <class Resource, class Identifier>
void ResourceHolder<Resource, Identifier>::load (Identifier id, std::unique_ptr<Resource> &resource) {
    auto inserted = resourceMap.insert(std::make_pair(id, std::move(resource)));
    assert(inserted.second);
}

template <class Resource, class Identifier>
void ResourceHolder<Resource, Identifier>::load (Identifier id, const std::string& filename) {
    std::unique_ptr<Resource> resource(std::make_unique<Resource>());
    if (!resource->loadFromFile(filename)) {
        throw std::runtime_error("ResourceHolder::load - Failed to load " + filename);
    }
    load(id, resource);
}

template <class Resource, class Identifier>
template <class Parameter>
void ResourceHolder<Resource, Identifier>::load (Identifier id, const std::string& filename, const Parameter& secondParam) {
    std::unique_ptr<Resource> resource(std::make_unique<Resource>());
    if (!resource->loadFromFile(filename, secondParam)) {
        throw std::runtime_error("ResourceHolder::load - Failed to load " + filename);
    }
    load(id, resource);
}

template <class Resource, class Identifier>
Resource &ResourceHolder<Resource, Identifier>::get (Identifier id) {
    auto found = resourceMap.find(id);
    assert(found != resourceMap.end());
    return *(found->second);
}

template <class Resource, class Identifier>
const Resource &ResourceHolder<Resource, Identifier>::get (Identifier id) const {
    auto found = resourceMap.find(id);
    assert(found != resourceMap.end());
    return *(found->second);
}

template <class Resource, class Identifier>
Resource &ResourceHolder<Resource, Identifier>::operator[] (Identifier id) {
    return get(id);
}

template <class Resource, class Identifier>
const Resource &ResourceHolder<Resource, Identifier>::operator[] (Identifier id) const {
    return get(id);
}

#endif